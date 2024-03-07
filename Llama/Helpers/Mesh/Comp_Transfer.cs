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


namespace Llama.Helpers.Mesh
{
    /// <summary>
    /// A grasshopper component trimming a mesh with boundary polylines.
    /// </summary>
    public class Comp_Transfer : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_Transfer"/> class.
        /// </summary>
        public Comp_Transfer()
          : base("Transfer", "Transfer",
              "Transfer points from target's flat configuration to its initial configuration.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points on the flat configuration to transfer", GH_Kernel.GH_ParamAccess.tree);

            pManager.AddMeshParameter("Flat Target", "Tf", "Target mesh in its flat configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddMeshParameter("Initial Target", "Ti", "Target mesh in its initial configuration.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Transfered Points", "T", "Points transfered on the initial configuration.", GH_Kernel.GH_ParamAccess.tree);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Point> struct_Points = new GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Point>();

            RH_Geo.Mesh targetFlat = new RH_Geo.Mesh();
            RH_Geo.Mesh targetInitial = new RH_Geo.Mesh();

            // ----- Get Inputs ----- //

            if (!DA.GetDataTree(0, out struct_Points)) { return; } ;

            if (!DA.GetData(1, ref targetFlat)) { return; } ;
            if (!DA.GetData(2, ref targetInitial)) { return; } ;

            // ----- Core ----- //

            GH.DataTree<RH_Geo.Point3d> points = new GH.DataTree<RH_Geo.Point3d>();

            for (int i = 0; i < struct_Points.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = struct_Points.Paths[i];
                List<GH_Kernel.Types.GH_Point> branch = struct_Points[path];

                for (int j = 0; j < branch.Count; j++)
                {
                    RH_Geo.Point3d position = branch[j].Value;
                    RH_Geo.MeshPoint point = targetFlat.ClosestMeshPoint(position, 0);

                    RH_Geo.MeshFace face = targetInitial.Faces[point.FaceIndex];
                    if (!face.IsTriangle)
                    {
                        this.AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Error, "The barycentric coordinate method is not implemented for none triangular faces.");
                        return;
                    }

                    RH_Geo.Point3d a = targetInitial.Vertices[face.A];
                    RH_Geo.Point3d b = targetInitial.Vertices[face.B];
                    RH_Geo.Point3d c = targetInitial.Vertices[face.C];

                    double weightA = point.T[0];
                    double weightB = point.T[1];
                    double weightC = point.T[2];

                    RH_Geo.Point3d transfered = (weightA * a) + (weightB * b) + (weightC * c);

                    points.Add(transfered, path);
                }
            }

            // ----- Set Output ----- //

            DA.SetDataTree(0, points);

        }

        #endregion  

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{29019D72-B4C5-4257-BE49-D78B1A831095}");

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
