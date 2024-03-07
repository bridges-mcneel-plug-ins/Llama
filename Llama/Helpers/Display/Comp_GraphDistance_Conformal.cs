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
    public class Comp_GraphDistance_Conformal : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_GraphDistance_Conformal"/> class.
        /// </summary>
        public Comp_GraphDistance_Conformal()
          : base("Graph Distance 2", "G. Dist. 2",
              "Transfer a mesh and Sort the edges of both mesh according to their graph distance to specified points.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Pf", "Pattern mesh whose edges needs to be sorted.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddMeshParameter("Mesh", "Tf", "Target mesh in its flat configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddMeshParameter("Mesh", "Ti", "Target mesh in its initial configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddPointParameter("Reference Points", "P", "Points to compute the graph distance", GH_Kernel.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Edges", "Ef", "Edges of the flat pattern sorted by graph distance.", GH_Kernel.GH_ParamAccess.tree);
            pManager.AddLineParameter("Edges", "Ep", "Edges of the projected pattern sorted by graph distance.", GH_Kernel.GH_ParamAccess.tree);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            RH_Geo.Mesh pattern = new RH_Geo.Mesh();
            RH_Geo.Mesh flat = new RH_Geo.Mesh();
            RH_Geo.Mesh initial = new RH_Geo.Mesh();
            List<RH_Geo.Point3d> points = new List<RH_Geo.Point3d>();

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref pattern)) { return; };
            if (!DA.GetData(1, ref flat)) { return; };
            if (!DA.GetData(2, ref initial)) { return; };
            if (!DA.GetDataList(3, points)) { return; };

            // ----- Transfer Mesh ----- //

            RH_Geo.PointCloud cloud = new RH_Geo.PointCloud();

            //RH_Geo.Mesh projected = new RH_Geo.Mesh();

            He.Mesh<Euc3D.Point> hePattern = new He.Mesh<Euc3D.Point>();
            He.Mesh<Euc3D.Point> heProjected = new He.Mesh<Euc3D.Point>();


            Dictionary<int, int> i_OldToNew = new Dictionary<int, int>(pattern.Vertices.Count);

            foreach(RH_Geo.MeshNgon ngon in pattern.GetNgonAndFacesEnumerable())
            {
                uint[] u_FaceVertices = ngon.BoundaryVertexIndexList();

                // Remove duplicate points in (true) ngon faces
                if (u_FaceVertices[0] == u_FaceVertices[u_FaceVertices.Length - 1])
                {
                    uint[] u_Temporary = new uint[u_FaceVertices.Length - 1];
                    for (int i = 0; i < u_Temporary.Length; i++) { u_Temporary[i] = u_FaceVertices[i]; }

                    u_FaceVertices = u_Temporary;
                }

                // Create the list of face vertices
                List<int> i_FaceVertices = new List<int>(u_FaceVertices.Length);
                for (int i = 0; i < u_FaceVertices.Length; i++)
                {
                    int i_FaceTopoVertex = pattern.TopologyVertices.TopologyVertexIndex((int)u_FaceVertices[i]);

                    int i_HeVertex;
                    if (i_OldToNew.ContainsKey(i_FaceTopoVertex)) { i_HeVertex = i_OldToNew[i_FaceTopoVertex]; }
                    else 
                    {
                        // For hePattern
                        RH_Geo.Point3d point3d = (RH_Geo.Point3d)pattern.TopologyVertices[i_FaceTopoVertex];
                        point3d.CastTo(out Euc3D.Point point);

                        // For heProjected
                        RH_Geo.MeshPoint meshPoint = flat.ClosestMeshPoint(point3d, 1e-4);

                        RH_Geo.MeshFace face = initial.Faces[meshPoint.FaceIndex];
                        if (!face.IsTriangle) { throw new Exception("The barycentric coordinate method is not implemented for none triangular faces."); }

                        RH_Geo.Point3d a = initial.Vertices[face.A];
                        RH_Geo.Point3d b = initial.Vertices[face.B];
                        RH_Geo.Point3d c = initial.Vertices[face.C];

                        double weightA = meshPoint.T[0];
                        double weightB = meshPoint.T[1];
                        double weightC = meshPoint.T[2];

                        RH_Geo.Point3d transferedPoint3d = (weightA * a) + (weightB * b) + (weightC * c);
                        transferedPoint3d.CastTo(out Euc3D.Point transferedPoint);

                        // Add the vertex
                        cloud.Add(point3d);
                        He.Vertex<Euc3D.Point> vertex = hePattern.AddVertex(point);
                        He.Vertex<Euc3D.Point> transferedVertex = heProjected.AddVertex(transferedPoint);

                        // Manages helpers
                        i_OldToNew.Add(i_FaceTopoVertex, vertex.Index);
                        i_HeVertex = vertex.Index;
                    }

                    // Add Face
                    i_FaceVertices.Add(i_HeVertex);
                }

                hePattern.AddFace(i_FaceVertices);
                heProjected.AddFace(i_FaceVertices);
            }




            // ----- Get Closest Vertices ----- //

            List<He.Vertex<Euc3D.Point>> vertices = new List<He.Vertex<Euc3D.Point>>();
            for (int i = 0; i < points.Count; i++)
            {
                int i_Vertex = cloud.ClosestPoint(points[i]);
                vertices.Add(hePattern.GetVertex(i_Vertex));
            }

            // ----- Core ----- //

            bool[] isEdgeUsed = new bool[hePattern.EdgeCount];
            GH.DataTree<RH_Geo.Line> linesPattern = new GH.DataTree<RH_Geo.Line>();
            GH.DataTree<RH_Geo.Line> linesProjected = new GH.DataTree<RH_Geo.Line>();

            while (vertices.Count != 0)
            {
                GH_Kernel.Data.GH_Path path = new GH_Kernel.Data.GH_Path(linesPattern.BranchCount);
                List<He.Vertex<Euc3D.Point>> nextVertices = new List<He.Vertex<Euc3D.Point>>();

                for (int i = 0; i < vertices.Count; i++)
                {
                    IReadOnlyList<He.Halfedge<Euc3D.Point>> outgoingHalfedges = vertices[i].OutgoingHalfedges();

                    for (int j = 0; j < outgoingHalfedges.Count; j++)
                    {
                        He.Halfedge<Euc3D.Point> outgoingHalfedge = outgoingHalfedges[j];

                        int i_Edge = Math.DivRem(outgoingHalfedge.Index, 2, out _);
                        if (isEdgeUsed[i_Edge] == true) { continue; }

                        // For hePattern
                        outgoingHalfedge.StartVertex.Position.CastTo(out RH_Geo.Point3d start);
                        outgoingHalfedge.EndVertex.Position.CastTo(out RH_Geo.Point3d end);

                        RH_Geo.Line edgeLine = new RH_Geo.Line(start, end);

                        // For heProjected
                        He.Halfedge<Euc3D.Point> outgoingHalfedgeProjected = heProjected.GetHalfedge(outgoingHalfedge.Index);

                        outgoingHalfedgeProjected.StartVertex.Position.CastTo(out RH_Geo.Point3d startProjected);
                        outgoingHalfedgeProjected.EndVertex.Position.CastTo(out RH_Geo.Point3d endProjected);

                        RH_Geo.Line edgeLineProjected = new RH_Geo.Line(startProjected, endProjected);

                        // Add lines
                        linesPattern.Add(edgeLine, path);
                        linesProjected.Add(edgeLineProjected, path);

                        // Manages helpers
                        isEdgeUsed[i_Edge] = true;
                        nextVertices.Add(outgoingHalfedge.EndVertex);
                    }
                }

                vertices = nextVertices;
            }


            // ----- Set Output ----- //

            DA.SetDataTree(0, linesPattern);
            DA.SetDataTree(1, linesProjected);


        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{939D58D3-FBC4-415E-8E7B-9505291F1B8C}");

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
