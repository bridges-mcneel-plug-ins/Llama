using System;
using System.Collections.Generic;

using RH_Geo = Rhino.Geometry;

using GH = Grasshopper;
using GH_Kernel = Grasshopper.Kernel;

using He = BRIDGES.DataStructures.PolyhedralMeshes.HalfedgeMesh;
using Euc3D = BRIDGES.Geometry.Euclidean3D;

using BRIDGES.McNeel.Rhino.Extensions.Geometry.Euclidean3D;
using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;


namespace Llama.Helpers.Display
{
    /// <summary>
    /// A grasshopper component computing the graph distance of edges to reference points.
    /// </summary>
    internal class Comp_EdgeDistance : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_EdgeDistance"/> class.
        /// </summary>
        public Comp_EdgeDistance()
          : base("Graph Distance", "G. Dist.",
              "Sort the edges of a mesh according to their graph distance to specified points.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Mesh whose edges needs to be sorted.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddPointParameter("Reference Points", "P", "Points to compute the graph distance", GH_Kernel.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Edges", "E", "Edges sorted by graph distance.", GH_Kernel.GH_ParamAccess.tree);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            RH_Geo.Mesh mesh = new RH_Geo.Mesh();
            List<RH_Geo.Point3d> points = new List<RH_Geo.Point3d>();

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref mesh)) { return; } ;
            if (!DA.GetDataList(1, points)) { return; };

            // ----- Transfer Mesh ----- //

            RH_Geo.PointCloud cloud = new RH_Geo.PointCloud();

            He.Mesh<Euc3D.Point> heMesh = new He.Mesh<Euc3D.Point>();
            Dictionary<int, int> i_OldToNew = new Dictionary<int, int>(mesh.Vertices.Count);

            foreach(RH_Geo.MeshNgon ngon in mesh.GetNgonAndFacesEnumerable())
            {
                uint[] u_FaceVertices = ngon.BoundaryVertexIndexList();

                if (u_FaceVertices[0] == u_FaceVertices[u_FaceVertices.Length - 1])
                {
                    uint[] u_Temporary = new uint[u_FaceVertices.Length - 1];
                    for (int i = 0; i < u_Temporary.Length; i++) { u_Temporary[i] = u_FaceVertices[i]; }

                    u_FaceVertices = u_Temporary;
                }


                List<int> i_FaceVertices = new List<int>(u_FaceVertices.Length);
                for (int i = 0; i < u_FaceVertices.Length; i++)
                {
                    int i_FaceTopoVertex = mesh.TopologyVertices.TopologyVertexIndex((int)u_FaceVertices[i]);

                    int i_HeVertex;
                    if (i_OldToNew.ContainsKey(i_FaceTopoVertex)) { i_HeVertex = i_OldToNew[i_FaceTopoVertex]; }
                    else 
                    {
                        RH_Geo.Point3d point3d = (RH_Geo.Point3d)mesh.TopologyVertices[i_FaceTopoVertex];
                        point3d.CastTo(out Euc3D.Point point);

                        cloud.Add(point3d);
                        He.Vertex<Euc3D.Point> vertex = heMesh.AddVertex(point);

                        i_OldToNew.Add(i_FaceTopoVertex, vertex.Index);
                        i_HeVertex = vertex.Index;
                    }
                    i_FaceVertices.Add(i_HeVertex);
                }

                heMesh.AddFace(i_FaceVertices);
            }

            // ----- Get Closest Vertices ----- //

            List<He.Vertex<Euc3D.Point>> vertices = new List<He.Vertex<Euc3D.Point>>();
            for (int i = 0; i < points.Count; i++)
            {
                int i_Vertex = cloud.ClosestPoint(points[i]);
                vertices.Add(heMesh.GetVertex(i_Vertex));
            }

            // ----- Core ----- //

            bool[] isEdgeUsed = new bool[heMesh.EdgeCount];
            GH.DataTree<RH_Geo.Line> lines = new GH.DataTree<RH_Geo.Line>();

            while (vertices.Count != 0)
            {
                GH_Kernel.Data.GH_Path path = new GH_Kernel.Data.GH_Path(lines.BranchCount);
                List<He.Vertex<Euc3D.Point>> nextVertices = new List<He.Vertex<Euc3D.Point>>();

                for (int i = 0; i < vertices.Count; i++)
                {
                    IReadOnlyList<He.Halfedge<Euc3D.Point>> outgoingHalfedges = vertices[i].OutgoingHalfedges();

                    for (int j = 0; j < outgoingHalfedges.Count; j++)
                    {
                        He.Halfedge<Euc3D.Point> outgoingHalfedge = outgoingHalfedges[j];

                        int i_Edge = Math.DivRem(outgoingHalfedge.Index, 2, out _);
                        if (isEdgeUsed[i_Edge] == true) { continue; }

                        outgoingHalfedge.StartVertex.Position.CastTo(out RH_Geo.Point3d start);
                        outgoingHalfedge.EndVertex.Position.CastTo(out RH_Geo.Point3d end);

                        RH_Geo.Line edgeLine = new RH_Geo.Line(start, end);

                        lines.Add(edgeLine, path);
                        isEdgeUsed[i_Edge] = true;

                        nextVertices.Add(outgoingHalfedge.EndVertex);
                    }
                }

                vertices = nextVertices;
            }


            // ----- Set Output ----- //

            DA.SetDataTree(0, lines);


        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{92AF524B-D602-499E-AAAA-627C0F81CB15}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return null;
            }
        }

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Display;


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }
}
