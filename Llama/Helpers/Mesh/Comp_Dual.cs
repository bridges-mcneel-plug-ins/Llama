using System;
using System.Collections.Generic;

using RH = Rhino;
using RH_Geo = Rhino.Geometry;

using GH = Grasshopper;
using GH_Kernel = Grasshopper.Kernel;

using He = BRIDGES.DataStructures.PolyhedralMeshes.HalfedgeMesh;
using Euc3D = BRIDGES.Geometry.Euclidean3D;

using BRIDGES.McNeel.Rhino.Extensions.Geometry.Euclidean3D;
using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;
using System.IO;


namespace Llama.Helpers.Mesh
{
    /// <summary>
    /// A grasshopper component trimming a mesh with boundary polylines.
    /// </summary>
    public class Comp_Dual : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_Dual"/> class.
        /// </summary>
        public Comp_Dual()
          : base("Dual Mesh", "Dual",
              "Create tje dual of a mesh.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Vertices", "V", "Vertices of the flat trimmed pattern.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddIntegerParameter("Faces", "F", "Faces of the flat trimmed pattern.", GH_Kernel.GH_ParamAccess.tree);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Vertices", "V", "Vertices of the dual.", GH_Kernel.GH_ParamAccess.tree);
            pManager.AddPathParameter("Faces Vertices", "F", "Faces vertices of the dual.", GH_Kernel.GH_ParamAccess.tree);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            List<RH_Geo.Point3d> vertices = new List<RH_Geo.Point3d>();
            GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Integer> gh_Structure = new GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Integer>();

            // ----- Get Inputs ----- //

            if (!DA.GetDataList(0, vertices)) { return; } ;
            if (!DA.GetDataTree(1, out gh_Structure)) { return; } ;

            #region GH_Structure to DataTree

            GH.DataTree<int> i_FacesInternalBorders = new GH.DataTree<int>();
            for (int i = 0; i < gh_Structure.PathCount; i++)
            {
                GH_Kernel.Data.GH_Path path = gh_Structure.Paths[i];
                List<GH_Kernel.Types.GH_Integer> branch = gh_Structure[path];
                for (int j = 0; j < branch.Count; j++)
                {
                    i_FacesInternalBorders.Add(branch[j].Value, path);
                }
            }

            #endregion

            // ----- Core ----- //

            #region Create Halfedge Mesh from Inputs

            He.Mesh<Euc3D.Point> hePrimal = new He.Mesh<Euc3D.Point>();
            Dictionary<GH_Kernel.Data.GH_Path, int> i_FaceToHeFace = new Dictionary<GH_Kernel.Data.GH_Path, int>(i_FacesInternalBorders.Paths.Count);
            List<GH_Kernel.Data.GH_Path> i_HeFaceToFace = new List<GH_Kernel.Data.GH_Path>(i_FacesInternalBorders.Paths.Count);

            for (int i = 0; i < vertices.Count; i++)
            {
                vertices[i].CastTo(out Euc3D.Point position);
                hePrimal.AddVertex(position);
            }
            for (int i = 0; i < i_FacesInternalBorders.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = i_FacesInternalBorders.Paths[i];

                List<int> i_FaceInternalBorder = i_FacesInternalBorders.Branch(path);
                He.Face<Euc3D.Point> face = hePrimal.AddFace(i_FaceInternalBorder);

                i_FaceToHeFace.Add(path, face.Index);
                i_HeFaceToFace.Add(path);
            }

            #endregion

            #region Create Vertices of the Dual

            List<RH_Geo.Point3d> dualVertices = new List<RH_Geo.Point3d>();
            GH.DataTree<RH_Geo.Point3d> dualVertices_Tree = new GH.DataTree<RH_Geo.Point3d>();

            for (int i = 0; i < i_FacesInternalBorders.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = i_FacesInternalBorders.Paths[i];

                List<int> i_FaceInternalBorder = i_FacesInternalBorders.Branch(path);

                RH_Geo.Point3d centroid = new RH_Geo.Point3d(0d, 0d, 0d);
                for (int j = 0; j < i_FaceInternalBorder.Count; j++)
                {
                    RH_Geo.Point3d position = vertices[i_FaceInternalBorder[j]];
                    centroid += position;
                }
                centroid /= i_FaceInternalBorder.Count;

                dualVertices.Add(centroid);
                dualVertices_Tree.Add(centroid, path);
            }

            #endregion

            #region Create Faces of the Dual

            GH.DataTree<int> i_DualFacesVertices = new GH.DataTree<int>();
            GH.DataTree<GH_Kernel.Data.GH_Path> p_DualFacesVertices = new GH.DataTree<GH_Kernel.Data.GH_Path>();
            for (int i = 0; i < hePrimal.VertexCount; i++)
            {
                He.Vertex<Euc3D.Point> heVertex = hePrimal.GetVertex(i);
                if(heVertex.IsBoundary()) { continue; }

                var adjacentFaces = heVertex.AdjacentFaces();

                List<int> i_DualFaceVertices = new List<int>();
                List<GH_Kernel.Data.GH_Path> p_DualFaceVertices = new List<GH_Kernel.Data.GH_Path>();
                for (int j = 0; j < adjacentFaces.Count; j++)
                {
                    int i_HeFace = adjacentFaces[j].Index;
                    i_DualFaceVertices.Add(i_HeFace);
                    p_DualFaceVertices.Add(i_HeFaceToFace[i_HeFace]);
                }

                i_DualFacesVertices.AddRange(i_DualFaceVertices, new GH_Kernel.Data.GH_Path(i));
                p_DualFacesVertices.AddRange(p_DualFaceVertices, new GH_Kernel.Data.GH_Path(i));
            }

            #endregion

            // ----- Set Output ----- //

            DA.SetDataTree(0, dualVertices_Tree);
            DA.SetDataTree(1, p_DualFacesVertices);

        }

        #endregion  

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{7EFD86B9-B695-4BD0-9DD9-2717A0D2184C}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return null;
            }
        }

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Mesh;


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

}
