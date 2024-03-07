using System;
using System.Linq;
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
    public class TrimAndJoin : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="TrimAndJoin"/> class.
        /// </summary>
        public TrimAndJoin()
          : base("Trim and Join", "T & J",
              "Trim a mesh with a set of boundary polylines defining an connected domain, and join the boundary segments in the mesh.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Mesh to operate on.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddCurveParameter("Boundary Polylines", "C", "Boundary polylines defining an connected domain.", GH_Kernel.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Trimmed Mesh", "T", "Mesh trimmed with the specified boundary polylines.", GH_Kernel.GH_ParamAccess.item);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            double tolerance = 1e-4;

            RH_Geo.Mesh mesh = new RH_Geo.Mesh();
            List<RH_Geo.Curve> curves = new List<RH_Geo.Curve>();

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref mesh)) { return; } ;
            if (!DA.GetDataList(1, curves)) { return; };


            // ----- Initialisation ----- //

            RH_Geo.MeshNgon[] ngonAndFaces = mesh.GetNgonAndFacesEnumerable().ToArray();


            // ----- Prepare Connected Domain ----- //

            RH_Geo.PolylineCurve[] boundaries = new RH_Geo.PolylineCurve[curves.Count];
            for (int i = 0; i < boundaries.Length; i++) 
            {
                if (curves[i].TryGetPolyline(out RH_Geo.Polyline polyline)) { boundaries[i] = polyline.ToPolylineCurve(); }
                else { throw new ArgumentException("The boundary curves must be polylines.");  }

            }

            RH_Geo.Brep[] breps = RH_Geo.Brep.CreatePlanarBreps(boundaries, tolerance);
            if (breps.Length != 1) { throw new ArgumentException("The boundary curves must define only one connected domain."); }
            RH_Geo.Brep domain = breps[0];

            RH_Geo.Vector3d domainNormal = domain.Faces[0].NormalAt(0.5, 0.5);

            // ----- Split Edges ----- //

            List<RH_Geo.Point3d> vertices;
            
            GH.DataTree<(int, double)> bordersSplitInfos;

            GH.DataTree<(int, int)> i_EdgesInternalSegments = EdgesTrimming(mesh, domain, boundaries, tolerance, out vertices, out bordersSplitInfos);

            // ----- Split Borders ----- //

            Dictionary<(int, int), RH_Geo.Polyline> boundariesSegments = BordersSplitting(boundaries, bordersSplitInfos);

            // ----- Face Internal Borders ----- //

            GH.DataTree<int> i_FacesInternalBorders = FacesInternalBorders(i_EdgesInternalSegments, mesh, ngonAndFaces);




            // ----- Face Internal Borders ----- //

            RH_Geo.Mesh trimmedMesh = new RH_Geo.Mesh();

            //----- Add Vertices -----/

            for (int i_V = 0; i_V < vertices.Count; i_V++)
            {
                trimmedMesh.Vertices.Add(vertices[i_V]);
            }

            // ----- Get faces Borders ----- //
            GH.DataTree<int> i_FacesBorders = FacesBorders(boundariesSegments, i_FacesInternalBorders, ref trimmedMesh, mesh, ngonAndFaces, tolerance);

            // ----- Add Faces ----- //

            for (int i = 0; i < i_FacesBorders.BranchCount; i++)
            {
                List<int> i_FaceBorder = i_FacesBorders.Branch(i);
                AddFace(i_FaceBorder, ref trimmedMesh, domainNormal, tolerance);
            }

            // ----- Set Output ----- //

            DA.SetData(0, trimmedMesh);

        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{D777C772-323B-4418-84F9-81DCD9F0270B}");

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


        #region Edges Trimming

        /// <summary>
        /// Trims the edges of the mesh with the boundary polylines.
        /// </summary>
        /// <param name="mesh"> Mesh to operate on. </param>
        /// <param name="domain"> Brep representing the connected domain. </param>
        /// <param name="boundaries"> Curves defining the boundaries of admissible domain. </param>
        /// <param name="tolerance"> Precision of the computation. </param>
        /// <param name="vertices"> List of vertices for the new trimmed mesh. </param>
        /// <param name="boundariesSplitInfos"> Information to split the boundaries. (Index of the splitting vertex, Splitting parameter of the border). </param>
        /// <returns> The edges internal segments, usefull for the new trimmed mesh, defined by the indices of the end vertices. </returns>
        private static GH.DataTree<(int, int)> EdgesTrimming(RH_Geo.Mesh mesh, RH_Geo.Brep domain, RH_Geo.Curve[] boundaries, double tolerance,
          out List<RH_Geo.Point3d> vertices, out GH.DataTree<(int, double)> boundariesSplitInfos)
        {
            GH.DataTree<(int, int)> i_EdgesSegments = new GH.DataTree<ValueTuple<int, int>>();

            // Initialise the necessary datas
            vertices = new List<RH_Geo.Point3d>();
            boundariesSplitInfos = new GH.DataTree<ValueTuple<int, double>>();

            Dictionary<int, int> i_OldToNew = new Dictionary<int, int>();

            // Iterate on the edges of the pattern

            int edgeCount = mesh.TopologyEdges.Count;
            for (int i = 0; i < edgeCount; i++)
            {
                RH_Geo.Line edgeLine = mesh.TopologyEdges.EdgeLine(i);
                RH_Geo.LineCurve edgeCurve = new RH_Geo.LineCurve(edgeLine);

                RH.IndexPair i_EndVertices = mesh.TopologyEdges.GetTopologyVertices(i);

                // ----- Compute Intersections ----- //

                List<(int, double, double)> edgeSplitParameters = EdgeSplitParameters(edgeCurve, boundaries, tolerance);

                // ----- Translate Points to Vertex Indices ----- //
                
                List<int> i_EdgeVertices = new List<int>();


                RH_Geo.Point3d projected = domain.ClosestPoint(edgeLine.From);
                bool isStartInside = edgeLine.From.DistanceTo(projected) < tolerance;

                // Manage Start Vertex
                if (isStartInside)
                {
                    int newIndex = ManageExistingVertex(i_EndVertices.I, mesh, ref vertices, ref i_OldToNew);
                    i_EdgeVertices.Add(newIndex);
                }

                // Manage Intersection Vertices
                foreach ((int, double, double) edgeSplitParameter in edgeSplitParameters)
                {
                    int i_Vertex = vertices.Count;
                    i_EdgeVertices.Add(i_Vertex);

                    RH_Geo.Point3d vertex = edgeCurve.PointAt(edgeSplitParameter.Item3);
                    vertices.Add(vertex);

                    boundariesSplitInfos.Add((i_Vertex, edgeSplitParameter.Item2), new GH_Kernel.Data.GH_Path(edgeSplitParameter.Item1));
                }

                // Manage End Vertex
                if ((isStartInside & edgeSplitParameters.Count % 2 == 0) | (!isStartInside & edgeSplitParameters.Count % 2 == 1))
                {
                    int newIndex = ManageExistingVertex(i_EndVertices.J, mesh, ref vertices, ref i_OldToNew);
                    i_EdgeVertices.Add(newIndex);
                }

                // ----- Sort Splitted|Internal|External Edges ----- //

                if (i_EdgeVertices.Count == 0) { i_EdgesSegments.AddRange(new List<ValueTuple<int, int>>(), new GH_Kernel.Data.GH_Path(i)); }
                else
                {
                    for (int i_EV = 1; i_EV < i_EdgeVertices.Count; i_EV += 2)
                    {
                        i_EdgesSegments.Add(ValueTuple.Create(i_EdgeVertices[i_EV - 1], i_EdgeVertices[i_EV]), new GH_Kernel.Data.GH_Path(i));
                    }
                }
                
                

                /*
                                for (int j = 1; j < i_EdgeVertices.Count; j += 2)
                                {
                                    i_EdgesSegments.Add((i_EdgeVertices[j - 1], i_EdgeVertices[j]), new GH_Kernel.Data.GH_Path(i));
                                }*/
            }

            return i_EdgesSegments;
        }


        /// <summary>
        /// Compute the intersection of the specified edge with the boundary curves.
        /// </summary>
        /// <param name="edgeCurve"> LineCurve representing the edge to intersect. </param>
        /// <param name="borders"> Curves to intersect the edge with. </param>
        /// <param name="tolerance"> Precision of the computation. </param>
        /// <returns>
        /// A list containing information for each intersection (BoundaryIndex, BoundaryParameter, EdgeParameter), sorted along the edge.
        /// </returns>
        private static List<(int, double, double)> EdgeSplitParameters(RH_Geo.LineCurve edgeCurve, RH_Geo.Curve[] borders, double tolerance)
        {
            List<ValueTuple<int, double, double>> edgeSplitParameters = new List<ValueTuple<int, double, double>>();

            for (int i_Bo = 0; i_Bo < borders.Length; i_Bo++)
            {
                RH_Geo.Intersect.CurveIntersections curveIntersections = RH_Geo.Intersect.Intersection.CurveCurve(borders[i_Bo], edgeCurve, tolerance, 1e-2 * tolerance);
                foreach (RH_Geo.Intersect.IntersectionEvent intersection in curveIntersections)
                {
                    if (intersection.IsOverlap) { throw new NotImplementedException("Edge-Border overlapping intersections are not managed."); }
                    if (intersection.IsPoint)
                    {
                        edgeSplitParameters.Add(ValueTuple.Create(i_Bo, intersection.ParameterA, intersection.ParameterB));
                    }
                }
            }

            edgeSplitParameters.Sort((x, y) => x.Item3.CompareTo(y.Item3));

            return edgeSplitParameters;
        }

        /// <summary>
        /// Gets the index of an existing vertex in the new trimmmed mesh. Adds it to the new trimmmed mesh if it was not already in.
        /// </summary>
        /// <param name="i_Vertex"> Index of the existing vertex in the initial mesh. </param>
        /// <param name="mesh"> Initial mesh. </param>
        /// <param name="newVertices"> List of new vertices og the new trimmed mesh. </param>
        /// <param name="i_OldToNew"> Dictionary defining the conversion of the vertex indices from the initial mesh to the new trimmed mesh. </param>
        /// <returns> The index of the vertex in the new trimmmed mesh. </returns>
        private static int ManageExistingVertex(int i_Vertex, RH_Geo.Mesh mesh, ref List<RH_Geo.Point3d> newVertices,
          ref Dictionary<int, int> i_OldToNew)
        {
            int i_NewVertex;
            if (i_OldToNew.ContainsKey(i_Vertex))
            {
                i_NewVertex = i_OldToNew[i_Vertex];
            }
            else // Adds the old vertex to the new list
            {
                i_NewVertex = newVertices.Count;

                i_OldToNew.Add(i_Vertex, i_NewVertex);
                newVertices.Add(mesh.TopologyVertices[i_Vertex]);
            }

            return i_NewVertex;
        }

        #endregion

        #region Border Splitting

        /// <summary>
        /// Splits the boundaries of the domain according to its intersections with the pattern edges.
        /// </summary>
        /// <param name="boundaries"> Curves defining the boundaries of admissible domain. </param>
        /// <param name="boundariesSplitInfos"> Information to split the boundaries. (Index of the splitting vertex, Splitting parameter of the border). </param>
        /// <returns> The segments of the boundaries, identified by the indices of their end vertices in the new trimmed mesh. </returns>
        private static Dictionary<ValueTuple<int, int>, RH_Geo.Polyline> BordersSplitting(RH_Geo.PolylineCurve[] boundaries, GH.DataTree<ValueTuple<int, double>> boundariesSplitInfos)
        {
            // Initialise the necessary datas
            Dictionary<ValueTuple<int, int>, RH_Geo.Polyline> boundariesSegments = new Dictionary<ValueTuple<int, int>, RH_Geo.Polyline>();

            // Iteration
            for (int i = 0; i < boundaries.Length; i++)
            {
                RH_Geo.PolylineCurve boundary = boundaries[i];

                List<(int, double)> boundarySplitInfos = boundariesSplitInfos.Branch(i);
                boundarySplitInfos.Sort((x, y) => x.Item2.CompareTo(y.Item2));

                //----- Split the border -----//

                List<double> boundarySplitParameters = new List<double>(boundarySplitInfos.Count);
                for (int j = 0; j < boundarySplitInfos.Count; j++)
                {
                    boundarySplitParameters.Add(boundarySplitInfos[j].Item2);
                }

                RH_Geo.Curve[] boundarySegments = boundary.Split(boundarySplitParameters);

                //----- Stores the segments -----//

                for (int j = 0; j < boundarySegments.Length; j++)
                {
                    RH_Geo.Curve boundarySegment = boundarySegments[j];

                    RH_Geo.Polyline polyline;
                    if (boundarySegment.TryGetPolyline(out polyline))
                    {
                        int i_StartVertex = boundarySplitInfos[j].Item1;
                        int i_EndVertex = boundarySplitInfos[(j + 1) % boundarySplitInfos.Count].Item1;

                        boundariesSegments.Add((i_StartVertex, i_EndVertex), polyline);
                    }
                    else { throw new ArgumentException("The boundary curves must be polylines so that a trimmed mesh can be rebuilt."); }
                }
            }

            return boundariesSegments;
        }

        #endregion

        #region Faces Internal Borders

        /// <summary>
        /// Joins the edge segments of each initial face, and create the border of the trimmed face(s).
        /// </summary>
        /// <param name="i_EdgesInternalSegments"> Edges internal segments usefull for the new trimmed mesh, defined by the indices of the end vertices. </param>
        /// <param name="mesh"> Mesh being trimmed. </param>
        /// <param name="ngonAndFaces"> Ngon and face of the mesh being trimmed. </param>
        /// <returns> Indices of the face vertices, forming the internal border of each face. </returns>
        private static GH.DataTree<int> FacesInternalBorders(GH.DataTree<(int, int)> i_EdgesInternalSegments, RH_Geo.Mesh mesh, RH_Geo.MeshNgon[] ngonAndFaces)
        {
            // Initialise the necessary datas
            GH.DataTree<int> i_FacesInternalBorders = new GH.DataTree<int>();

            // Iterate on each faces

            for (int i = 0; i < ngonAndFaces.Length; i++)
            {
                RH_Geo.MeshNgon ngon = ngonAndFaces[i];

                //----- Get the ngon edges -----//
                uint[] u_FaceVertices = ngon.BoundaryVertexIndexList();
                int faceVertexCount = u_FaceVertices.Length;

                int[] i_FaceEdges = new int[faceVertexCount];
                for (int j = 0; j < faceVertexCount; j++)
                {
                    int i_Start = Convert.ToInt32(u_FaceVertices[j]);
                    int i_End = Convert.ToInt32(u_FaceVertices[(j + 1) % faceVertexCount]);

                    int i_TopoStart = mesh.TopologyVertices.TopologyVertexIndex(i_Start);
                    int i_TopoEnd = mesh.TopologyVertices.TopologyVertexIndex(i_End);

                    i_FaceEdges[j] = mesh.TopologyEdges.GetEdgeIndex(i_TopoStart, i_TopoEnd);
                }

                //----- Evaluate whether the edge is flipped. -----//

                bool[] isFaceEdgeFlipped = FaceEdgesOrientation(i_FaceEdges, mesh);

                //----- Combine Face Segments -----//

                List<(int, int)> i_FaceInternalSegments = FaceInternalSegments(i_FaceEdges, isFaceEdgeFlipped, i_EdgesInternalSegments);

                //----- Create Face Borders -----//

                if (i_FaceInternalSegments.Count == 0) { continue; }

                GH.DataTree<int> i_FaceInternalBorders = FaceInternalBorders(i_FaceInternalSegments);

                //----- Stors the Face Borders -----//

                for (int j = 0; j < i_FaceInternalBorders.BranchCount; j++)
                {
                    i_FacesInternalBorders.AddRange(i_FaceInternalBorders.Branch(j), new GH_Kernel.Data.GH_Path(new int[2] { i, j }));
                }
            }

            return i_FacesInternalBorders;
        }


        /// <summary>
        /// Evaluates which face edges needs to be flipped in order for them to have a similar orientation.
        /// </summary>
        /// <param name="i_FaceEdges"> Indices of the face edges. </param>
        /// <param name="mesh"> Mesh to trim. </param>
        /// <returns> A boolean array evaluating which face edges needs to be flipped. </returns>
        private static bool[] FaceEdgesOrientation(int[] i_FaceEdges, RH_Geo.Mesh mesh)
        {
            bool[] isFlipped = new bool[i_FaceEdges.Length];

            // Initialise
            RH.IndexPair i_FirstEndVertices = mesh.TopologyEdges.GetTopologyVertices(i_FaceEdges[0]);
            RH.IndexPair i_SecondEndVertices = mesh.TopologyEdges.GetTopologyVertices(i_FaceEdges[1]);

            if (i_FirstEndVertices.J == i_SecondEndVertices.I) { }
            else if (i_FirstEndVertices.J == i_SecondEndVertices.J) { isFlipped[1] = true; }
            else if (i_FirstEndVertices.I == i_SecondEndVertices.I) { isFlipped[0] = true; }
            else { isFlipped[0] = true; isFlipped[1] = true; }

            // Iterate
            for (int i_FE = 2; i_FE < i_FaceEdges.Length; i_FE++)
            {
                RH.IndexPair prev = mesh.TopologyEdges.GetTopologyVertices(i_FaceEdges[i_FE - 1]);
                RH.IndexPair curr = mesh.TopologyEdges.GetTopologyVertices(i_FaceEdges[i_FE]);
                if ((!isFlipped[i_FE - 1] & prev.J != curr.I) | (isFlipped[i_FE - 1] & prev.I != curr.I))
                {
                    isFlipped[i_FE] = true;
                }
            }

            return isFlipped;
        }

        /// <summary>
        /// Combines the consecutive internal segments of the face.
        /// </summary>
        /// <param name="i_FaceEdges"> Indices of the face edges. </param>
        /// <param name="isFaceEdgeFlipped"> Boolean array evaluating which face internal edges needs to be flipped.</param>
        /// <param name="i_EdgesInternalSegments"> Edges segments usefull for the new trimmed mesh, defined by the indices of the end vertices. </param>
        /// <returns> A list of the face internal segments, defined by the indices of the end vertices. </returns>
        private static List<(int, int)> FaceInternalSegments(int[] i_FaceEdges, bool[] isFaceEdgeFlipped, GH.DataTree<(int, int)> i_EdgesInternalSegments)
        {
            List<(int, int)> i_FaceInternalSegments = new List<ValueTuple<int, int>>();

            for (int i_FE = 0; i_FE < i_FaceEdges.Length; i_FE++)
            {
                List<(int, int)> edgeSegments = new List<(int, int)>(i_EdgesInternalSegments.Branch(i_FaceEdges[i_FE]));

                if (isFaceEdgeFlipped[i_FE])
                {
                    edgeSegments.Reverse();
                    for (int i_ES = 0; i_ES < edgeSegments.Count; i_ES++)
                    {
                        (int, int) edgeSegment = edgeSegments[i_ES];
                        i_FaceInternalSegments.Add((edgeSegment.Item2, edgeSegment.Item1));
                    }
                }
                else { i_FaceInternalSegments.AddRange(edgeSegments); }
            }

            return i_FaceInternalSegments;
        }

        /// <summary>
        /// Joins the consecutive face internal segments to form the face internal borders.
        /// </summary>
        /// <param name="i_FaceInternalSegments"> List of the face internal segments, defined by the indices of the end vertices.</param>
        /// <returns> A list of the face internal borders, defined by the indices of its vertices. </returns>
        private static GH.DataTree<int> FaceInternalBorders(List<(int, int)> i_FaceInternalSegments)
        {
            GH.DataTree<int> i_FaceInternalBorders = new GH.DataTree<int>();

            // Initialise
            List<int> i_FaceInternalBorder = new List<int> { i_FaceInternalSegments[0].Item1, i_FaceInternalSegments[0].Item2 };

            // Iterate
            for (int i_FS = 1; i_FS < i_FaceInternalSegments.Count; i_FS++)
            {
                ValueTuple<int, int> i_FaceSegment = i_FaceInternalSegments[i_FS];
                if (i_FaceInternalBorder[i_FaceInternalBorder.Count - 1] == i_FaceSegment.Item1)
                {
                    i_FaceInternalBorder.Add(i_FaceSegment.Item2);
                }
                else
                {
                    i_FaceInternalBorders.AddRange(i_FaceInternalBorder, new GH_Kernel.Data.GH_Path(i_FaceInternalBorders.BranchCount));

                    i_FaceInternalBorder = new List<int> { i_FaceSegment.Item1, i_FaceSegment.Item2 };
                }
            }
            i_FaceInternalBorders.AddRange(i_FaceInternalBorder, new GH_Kernel.Data.GH_Path(i_FaceInternalBorders.BranchCount));

            // Finalise : Joins the last and the first borders if needed
            if (i_FaceInternalBorders.BranchCount > 1)
            {
                GH_Kernel.Data.GH_Path firstPath = new GH_Kernel.Data.GH_Path(0);
                GH_Kernel.Data.GH_Path lastPath = new GH_Kernel.Data.GH_Path(i_FaceInternalBorders.BranchCount - 1);

                List<int> first = i_FaceInternalBorders.Branch(firstPath);
                List<int> last = i_FaceInternalBorders.Branch(lastPath);

                if (first[0] == last[last.Count - 1])
                {
                    first.RemoveAt(0);
                    last.AddRange(first);
                    i_FaceInternalBorders.RemovePath(firstPath);
                }
            }

            return i_FaceInternalBorders;
        }

        #endregion

        #region Faces External Borders

        private static GH.DataTree<int> FacesBorders(Dictionary<(int, int), RH_Geo.Polyline> boundariesSegments, GH.DataTree<int> i_FacesInternalBorders, ref RH_Geo.Mesh trimmedMesh,
            RH_Geo.Mesh mesh, RH_Geo.MeshNgon[] ngonAndFaces, double tolerance)
        {
            GH.DataTree<int> i_FacesBorders = new GH.DataTree<int>();

            for (int i = 0; i < i_FacesInternalBorders.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = i_FacesInternalBorders.Paths[i];

                // The vertex's indices loops : the face is inside the admissible domain (and thus has no external border)
                List<int> i_FaceFirstInternalBorder = i_FacesInternalBorders.Branch(path);
                if (i_FaceFirstInternalBorder[0] == i_FaceFirstInternalBorder[i_FaceFirstInternalBorder.Count - 1])
                {
                    // ----- Store face border ----- //
                    for (int j = 0; j < i_FaceFirstInternalBorder.Count - 1; j++)
                    {
                        i_FacesBorders.Add(i_FaceFirstInternalBorder[j], path);
                    }
                }
                // Otherwise, the face has an external border
                else
                {
                    RH_Geo.Polyline intialFaceBorder = InitialFaceBorder(path, mesh, ngonAndFaces);

                    // ----- Retrieve the face borders ----- //
                    GH.DataTree<int> i_FaceInternalBorders = FaceInternalBorders(i, i_FacesInternalBorders);

                    // ----- Find the face external borders ----- //
                    GH.DataTree<(int, bool, int)> i_ExternalBorders = FaceExternalBorders(i_FaceInternalBorders, boundariesSegments, intialFaceBorder, tolerance);

                    // ----- Create face borders ----- //
                    GH.DataTree<int> i_FaceBorders = CreateFaceBorders(i_FaceInternalBorders, i_ExternalBorders, boundariesSegments, ref trimmedMesh);


                    // ----- Store face borders ----- //
                    for (int j = 0; j < i_FaceBorders.BranchCount; j++)
                    {
                        List<int> i_FaceBorder = i_FaceBorders.Branch(j);

                        int faceIndex = path.Indices[0];

                        i_FacesBorders.AddRange(i_FaceBorder, new GH_Kernel.Data.GH_Path(new int[2] { faceIndex, j }));
                    }

                    i += i_FaceInternalBorders.BranchCount - 1;
                }
            }

            return i_FacesBorders;
        }


        private static RH_Geo.Polyline InitialFaceBorder(GH_Kernel.Data.GH_Path path, RH_Geo.Mesh pattern, RH_Geo.MeshNgon[] ngonAndFaces)
        {
            // Get initial face border
            RH_Geo.MeshNgon ngon = ngonAndFaces[path.Indices[0]];
            uint[] u_FaceVertices = ngon.BoundaryVertexIndexList();

            int faceVertexCount = u_FaceVertices.Length;

            RH_Geo.Polyline faceTotalBorder = new RH_Geo.Polyline(faceVertexCount);
            for (int i_E = 0; i_E < faceVertexCount; i_E++)
            {
                int i_Start = Convert.ToInt32(u_FaceVertices[i_E]);
                int i_TopoStart = pattern.TopologyVertices.TopologyVertexIndex(i_Start);

                faceTotalBorder.Add(pattern.TopologyVertices[i_TopoStart]);
            }
            faceTotalBorder.Add(faceTotalBorder[0]);

            return faceTotalBorder;
        }


        /// <summary>
        /// Retrieves the internal borders of a face.
        /// </summary>
        /// <param name="i_P"> Path index of the face's first internal border. </param>
        /// <param name="i_FacesInternalBorders"> Indices of the face vertices, forming the internal border of each face. </param>
        /// <returns> A list of the face internal borders, defined by the indices of its vertices. </returns>
        private static GH.DataTree<int> FaceInternalBorders(int i_P, GH.DataTree<int> i_FacesInternalBorders)
        {
            GH_Kernel.Data.GH_Path path = i_FacesInternalBorders.Paths[i_P];

            // Get the path starting by "i_P"
            List<GH_Kernel.Data.GH_Path> paths = new List<GH_Kernel.Data.GH_Path> { path }; 

            while (i_P + 1 < i_FacesInternalBorders.Paths.Count && i_FacesInternalBorders.Paths[i_P + 1].Indices[0] == path.Indices[0])
            {
                paths.Add(i_FacesInternalBorders.Paths[i_P + 1]);
                i_P++;
            }

            // Get the face internal borders thanks to the selected paths
            GH.DataTree<int> i_FaceInternalBorders = new GH.DataTree<int>();
            for (int i = 0; i < paths.Count; i++)
            {
                i_FaceInternalBorders.AddRange(i_FacesInternalBorders.Branch(paths[i]), new GH_Kernel.Data.GH_Path(i));
            }

            return i_FaceInternalBorders;
        }
        

        private static GH.DataTree<(int, bool, int)> FaceExternalBorders(GH.DataTree<int> i_FaceInternalBorders, Dictionary<(int, int), RH_Geo.Polyline> boundariesSegments, 
            RH_Geo.Polyline initialfaceBorder, double tolerance)
        {
            RH_Geo.PolylineCurve crv_InitialfaceBorder = initialfaceBorder.ToPolylineCurve();

            // ----- Prepare ---- //

            int internalBordersEndsCount = 2 * i_FaceInternalBorders.BranchCount;

            bool[] isInternalBorderEndUsed = new bool[internalBordersEndsCount];
            int[] i_InternalBordersEnds = new int[internalBordersEndsCount];
            for (int i_B = 0; i_B < i_FaceInternalBorders.BranchCount; i_B++)
            {
                List<int> i_FaceInternalBorder = i_FaceInternalBorders.Branch(i_B);

                i_InternalBordersEnds[(2 * i_B)] = i_FaceInternalBorder[0];
                i_InternalBordersEnds[(2 * i_B) + 1] = i_FaceInternalBorder[i_FaceInternalBorder.Count - 1];
            }

            // ----- Search ---- //

            int branchCount = 0;
            GH.DataTree<(int, bool, int)> i_ExternalBorders = new GH.DataTree<(int, bool, int)>();
            for (int i_B = 0; i_B < internalBordersEndsCount / 2; i_B++)
            {
                List<(int, bool, int)> i_Border = new List<(int, bool, int)>();

                if (isInternalBorderEndUsed[2 * i_B]) { continue; }
                else
                {
                    FindConsecutiveBorders((2 * i_B) + 1, i_InternalBordersEnds, isInternalBorderEndUsed, ref i_Border,
                      boundariesSegments, crv_InitialfaceBorder, tolerance);
                }

                for (int i = 0; i < i_Border.Count; i++)
                {
                    i_ExternalBorders.Add(i_Border[i], new GH_Kernel.Data.GH_Path(branchCount));
                }
                branchCount++;
            }

            return i_ExternalBorders;
        }

        /// <summary>
        /// Recursive fonction finding the external borders of following an internal borders for a face.
        /// </summary>
        /// <param name="i"> Index of the internal border end whose external border needs to be found. </param>
        /// <param name="i_BordersEnds"> Indices of the end vertices of each internal face borders. </param>
        /// <param name="isBorderEndUsed"> Evalutates whether the internal border ends are used by an external border. </param>
        /// <param name="i_ExternalBorder">
        /// List to recursively fill with the information of the next external border found. 
        /// (Index of the starting <paramref name="i_BordersEnds"/>, Evaluate whether the external <paramref name="boundariesSegments"/> must be flipped, Index of the ending <paramref name="i_BordersEnds"/>). </param>
        /// <param name="boundariesSegments"> The segments of the boundaries, identified by the indices of their end vertices in the new trimmed mesh. </param>
        /// <param name="intialFaceBorder"> Border of the face whose borders are being managed. </param>
        /// <param name="tolerance"> Precision of the computation.  </param>
        private static void FindConsecutiveBorders(int i, int[] i_BordersEnds, bool[] isBorderEndUsed, ref List<(int, bool, int)> i_ExternalBorder,
            Dictionary<(int, int), RH_Geo.Polyline> boundariesSegments, RH_Geo.Curve intialFaceBorder, double tolerance)
        {
            for (int increment = 1; increment < i_BordersEnds.Length; increment++)
            {
                // Possible other end of the external border
                int j = (i + increment) % (i_BordersEnds.Length);
                if (isBorderEndUsed[j]) { continue; }

                //----- Searches for an External Border between the two points -----//

                // External border is found and oriented has expected
                if (boundariesSegments.ContainsKey(ValueTuple.Create(i_BordersEnds[i], i_BordersEnds[j])))
                {
                    i_ExternalBorder.Add(ValueTuple.Create(i, false, j));
                }
                // External border is found but is not oriented has expected
                else if (boundariesSegments.ContainsKey(ValueTuple.Create(i_BordersEnds[j], i_BordersEnds[i])))
                {
                    i_ExternalBorder.Add(ValueTuple.Create(i, true, j));
                }
                // No external border found
                else { continue; }


                //----- Ensures that the External Border found is inside the face -----//

                ValueTuple<int, bool, int> tuple = i_ExternalBorder[i_ExternalBorder.Count - 1];

                ValueTuple<int, int> i_ExternalBorderEnds;
                if (tuple.Item2) { i_ExternalBorderEnds = ValueTuple.Create(i_BordersEnds[tuple.Item3], i_BordersEnds[tuple.Item1]); }
                else { i_ExternalBorderEnds = ValueTuple.Create(i_BordersEnds[tuple.Item1], i_BordersEnds[tuple.Item3]); }

                RH_Geo.Polyline externalBorder = boundariesSegments[i_ExternalBorderEnds];
                RH_Geo.PolylineCurve cc = externalBorder.ToPolylineCurve();
                RH_Geo.Point3d centerPoint = cc.PointAtNormalizedLength(0.5); // externalBorder.CenterPoint();

                bool isCenterPointInside = IsPointInsideCurve(centerPoint, intialFaceBorder, tolerance);
                if (!isCenterPointInside) { i_ExternalBorder.RemoveAt(i_ExternalBorder.Count - 1); continue; }

                //----- Stores the External faces and moves forward -----//

                isBorderEndUsed[tuple.Item1] = true; isBorderEndUsed[tuple.Item3] = true;

                int _ = Math.DivRem(tuple.Item3, 2, out int remainderJ);
                int k = remainderJ == 0 ? tuple.Item3 + 1 : tuple.Item3 - 1;

                if (isBorderEndUsed[k]) { break; }
                else
                {
                    FindConsecutiveBorders(k, i_BordersEnds, isBorderEndUsed, ref i_ExternalBorder, boundariesSegments, intialFaceBorder, tolerance);
                    break;
                }
            }
        }


        #region To sperate

        private static GH.DataTree<int> CreateFaceBorders(GH.DataTree<int> i_FaceInternalBorders, GH.DataTree<ValueTuple<int, bool, int>> i_ExternalBorders,
            Dictionary<ValueTuple<int, int>, RH_Geo.Polyline> bordersSegments, ref RH_Geo.Mesh trimmedMesh)
        {
            GH.DataTree<int> faceBorders = new GH.DataTree<int>();
            for (int i_B = 0; i_B < i_ExternalBorders.BranchCount; i_B++)
            {
                GH_Kernel.Data.GH_Path path = new GH_Kernel.Data.GH_Path(i_B);
                List<ValueTuple<int, bool, int>> i_ExternalBorder = i_ExternalBorders.Branch(i_B);

                for (int i_EB = 0; i_EB < i_ExternalBorder.Count; i_EB++)
                {
                    int j_EB = (i_EB + 1) % i_ExternalBorder.Count;

                    ValueTuple<int, bool, int> tuple = i_ExternalBorder[i_EB];

                    // Add Internal Border

                    int remainderI = 0;
                    int quotientI = Math.DivRem(tuple.Item1, 2, out remainderI);

                    List<int> branchI = new List<int>(i_FaceInternalBorders.Branch(quotientI));
                    if (remainderI == 0) { branchI.Reverse(); }
                    for (int i = 0; i < branchI.Count; i++)
                    {
                        faceBorders.Add(branchI[i], path);
                    }

                    // Add Polyline

                    int remainderJ = 0;
                    int quotientJ = Math.DivRem(tuple.Item3, 2, out remainderJ);
                    List<int> branchJ = new List<int>(i_FaceInternalBorders.Branch(quotientJ));
                    if (remainderJ == 1) { branchJ.Reverse(); }

                    int i_V = branchI[branchI.Count - 1];
                    int j_V = branchJ[0];

                    RH_Geo.Polyline polyline;
                    if (tuple.Item2)
                    {
                        polyline = new RH_Geo.Polyline(bordersSegments[ValueTuple.Create(j_V, i_V)]);
                        polyline.Reverse();
                    }
                    else { polyline = bordersSegments[ValueTuple.Create(i_V, j_V)]; }


                    for (int i = 1; i < polyline.Count - 1; i++)
                    {
                        int index = trimmedMesh.Vertices.Add(polyline[i]);
                        faceBorders.Add(index, path);
                    }
                }
                //faceBorders.Add(faceBorders[path, 0]);
            }

            return faceBorders;
        }


        #endregion

        #endregion


        /// <summary>
        /// Evaluates whether a point is inside a closed curve.
        /// </summary>
        /// <param name="point"> Point to evaluate. </param>
        /// <param name="curve"> Closed curve defining an inner and outer domaine. </param>
        /// <param name="tolerance"> Precision of the computation. </param>
        /// <returns> <see langword="true"/>  if the point is inside the curve, <see langword="false"/> otherwise. </returns>
        private static bool IsPointInsideCurve(RH_Geo.Point3d point, RH_Geo.Curve curve, double tolerance)
        {
            if (curve.ClosestPoint(point, out _, 1e-2 * tolerance)) { return true; }

            // First Direction
            RH_Geo.Line line = new RH_Geo.Line(point, RH_Geo.Vector3d.XAxis, 1d);

            RH_Geo.Intersect.CurveIntersections curveIntersections =
             RH_Geo.Intersect.Intersection.CurveLine(curve, line, tolerance, 1e-2 * tolerance);

            List<double> parameters = new List<double>();

            bool exited = false;
            foreach (RH_Geo.Intersect.IntersectionEvent intersection in curveIntersections)
            {
                if (intersection.IsOverlap) { exited = true; break; }
                else if (intersection.IsPoint && 1e-2 * tolerance < intersection.ParameterB)
                {
                    parameters.Add(intersection.ParameterB);
                }
            }
            if (!exited) { return (parameters.Count % 2) == 1; }

            // Other direction
            line = new RH_Geo.Line(point, RH_Geo.Vector3d.XAxis + RH_Geo.Vector3d.YAxis, 1d);

            curveIntersections =
              RH_Geo.Intersect.Intersection.CurveLine(curve, line, tolerance, 1e-2 * tolerance);

            parameters = new List<double>();

            exited = false;
            foreach (RH_Geo.Intersect.IntersectionEvent intersection in curveIntersections)
            {
                if (intersection.IsOverlap) { exited = true; break; }
                else if (intersection.IsPoint && 1e-2 * tolerance < intersection.ParameterB)
                {
                    parameters.Add(intersection.ParameterB);
                }
            }
            if (!exited) { return (parameters.Count % 2) == 1; }

            throw new Exception("Ray tracing failed");
        }


        /// <summary>
        /// Add a face to the mesh from the vertex indices.
        /// </summary>
        /// <param name="mesh"> Mesh in which the face should be added. </param>
        /// <param name="i_FaceVertices"> Vertex indices of the face to add. </param>
        /// <param name="normal"> Vector for the orientation of the face. </param>
        /// <param name="tolerance"> Brep representing the admissible domain for the pattern. </param>
        void AddFace(List<int> i_FaceVertices, ref RH_Geo.Mesh mesh, RH_Geo.Vector3d normal, double tolerance)
        {
            GH.DataTree<int> i_FacesVertices = new GH.DataTree<int>();

            if (i_FaceVertices.Count == 3)
            {
                i_FacesVertices.AddRange(new int[3] { i_FaceVertices[0], i_FaceVertices[1], i_FaceVertices[2] }, new GH_Kernel.Data.GH_Path(i_FacesVertices.BranchCount));
                //mesh.Faces.AddFace(i_FaceVertices[0], i_FaceVertices[1], i_FaceVertices[2]);
            }
            else if (i_FaceVertices.Count == 4)
            {
                i_FacesVertices.AddRange(new int[4] { i_FaceVertices[0], i_FaceVertices[1], i_FaceVertices[2], i_FaceVertices[3] }, new GH_Kernel.Data.GH_Path(i_FacesVertices.BranchCount));
                //mesh.Faces.AddFace(i_FaceVertices[0], i_FaceVertices[1], i_FaceVertices[2], i_FaceVertices[3]);
            }
            else if (4 < i_FaceVertices.Count)
            {
                // Create face boundary
                RH_Geo.Polyline polyline = new RH_Geo.Polyline(i_FaceVertices.Count);
                for (int i_FV = 0; i_FV < i_FaceVertices.Count; i_FV++)
                {
                    RH_Geo.Point3d point = mesh.Vertices[i_FaceVertices[i_FV]];
                    polyline.Add(point);
                }
                polyline.Add(polyline[0]);

                RH_Geo.PolylineCurve border = polyline.ToPolylineCurve();

                // Start Delaunay Triangulation

                List<int> i_Faces = new List<int>();
                List<int> i_Vertices = new List<int>(i_FaceVertices);

                while (i_Vertices.Count > 2)
                {
                    List<int> i_NewVertices = new List<int>();
                    for (int i = 0; i < i_Vertices.Count; i += 2)
                    {
                        int j = i + 1; int k = i + 2;

                        if (j == i_Vertices.Count) { i_NewVertices.Add(i_Vertices[i]); }
                        else
                        {
                            if (k == i_Vertices.Count) { k = 0; }

                            RH_Geo.Point3d first = mesh.Vertices[i_Vertices[i]];
                            RH_Geo.Point3d last = mesh.Vertices[i_Vertices[k]];

                            bool isNotSplitted = true;
                            RH_Geo.LineCurve lineCurve = new RH_Geo.LineCurve(first, last);
                            RH_Geo.Intersect.CurveIntersections curveIntersections = RH_Geo.Intersect.Intersection.CurveCurve(border, lineCurve, tolerance, 1e-2 * tolerance);
                            foreach(RH_Geo.Intersect.IntersectionEvent intersection in curveIntersections)
                            {
                                if(intersection.IsOverlap)
                                {
                                    if (tolerance < intersection.PointA.DistanceTo(first) & tolerance < intersection.PointA.DistanceTo(last)
                                        & tolerance < intersection.PointB.DistanceTo(first) & tolerance < intersection.PointB.DistanceTo(last))
                                    {
                                        isNotSplitted = false; break;
                                    }
                                }
                                else if(intersection.IsPoint && (tolerance < intersection.PointA.DistanceToSquared(first) & tolerance < intersection.PointA.DistanceToSquared(last))) 
                                { 
                                    isNotSplitted = false; break; 
                                }
                            }

                            RH_Geo.Point3d center = (first + last) / 2d;
                            bool isInside = IsPointInsideCurve(center, border, tolerance);

                            if (isNotSplitted && isInside)
                            {
                                i_FacesVertices.AddRange(new int[3] { i_Vertices[i], i_Vertices[j], i_Vertices[k] }, new GH_Kernel.Data.GH_Path(i_FacesVertices.BranchCount));

                                //int i_Face = mesh.Faces.AddFace(i_Vertices[i], i_Vertices[j], i_Vertices[k]);
                                //i_Faces.Add(i_Face);
                                i_NewVertices.Add(i_Vertices[i]);
                            }
                            else { i_NewVertices.Add(i_Vertices[i]); i = i - 1; }
                        }
                    }

                    i_Vertices = new List<int>(i_NewVertices);
                }
            }

            // ----- Manages the Orientation ----- //

            List<int> i_FirstFaceVertices = i_FacesVertices.Branch(0);
            RH_Geo.Vector3d a = mesh.Vertices[i_FirstFaceVertices[1]] - mesh.Vertices[i_FirstFaceVertices[0]];
            RH_Geo.Vector3d b = mesh.Vertices[i_FirstFaceVertices[2]] - mesh.Vertices[i_FirstFaceVertices[0]];

            RH_Geo.Vector3d c = RH_Geo.Vector3d.CrossProduct(a, b);

            bool isWellOriented = c * normal > 0d;

            // ----- Add the Faces ----- //
            if (i_FacesVertices.BranchCount == 1)
            {
                List<int> branch = i_FacesVertices.Branch(0);
                if (!isWellOriented) { branch.Reverse(); }

                if (branch.Count == 3) { mesh.Faces.AddFace(branch[0], branch[1], branch[2]); }
                if (branch.Count == 4) { mesh.Faces.AddFace(branch[0], branch[1], branch[2], branch[3]); }
            }
            else
            {
                List<int> i_Faces = new List<int>(i_FacesVertices.BranchCount);
                for (int i = 0; i < i_FacesVertices.BranchCount; i++)
                {
                    List<int> branch = i_FacesVertices.Branch(i);
                    if (!isWellOriented) { branch.Reverse(); }

                    int i_Face = mesh.Faces.AddFace(branch[0], branch[1], branch[2]);
                    i_Faces.Add(i_Face);
                }

                if (!isWellOriented) { i_FaceVertices.Reverse(); }
                RH_Geo.MeshNgon ngon = RH_Geo.MeshNgon.Create(i_FaceVertices, i_Faces);
                mesh.Ngons.AddNgon(ngon);
            }
        }

    }

}
