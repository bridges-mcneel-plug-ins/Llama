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
    internal class Comp_MeshTrim_2 : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_MeshTrim_2"/> class.
        /// </summary>
        public Comp_MeshTrim_2()
          : base("Mesh Trim 2", "Trim 2",
              "Trim a mesh with a set of boundary polylines defining an connected domain.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Flat Pattern", "Pf", "Pattern mesh to trim in its flat configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddMeshParameter("Initial Target", "Ti", "Target mesh in its initial configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddMeshParameter("Flat Target", "Tf", "Target mesh in its flat configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddCurveParameter("Boundary Polylines", "C", "Boundary polylines defining an connected domain.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddSurfaceParameter("Height Surface", "S", "Surface defining the height of the dual mesh.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Primal Lines", "P", "Lines from the primal mesh.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddCurveParameter("Triangulation Lines", "T", "Lines connecting the primal and the dual.", GH_Kernel.GH_ParamAccess.tree);
            pManager.AddCurveParameter("Dual Lines", "p", "Lines from the dual mesh.", GH_Kernel.GH_ParamAccess.list);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            double tolerance = 1e-4;
            RH_Geo.Vector3d globalNormal = RH_Geo.Vector3d.ZAxis;

            RH_Geo.Mesh flatPattern = new RH_Geo.Mesh();
            RH_Geo.Mesh initialTarget = new RH_Geo.Mesh();
            RH_Geo.Mesh flatTarget = new RH_Geo.Mesh();

            List<RH_Geo.Curve> curves = new List<RH_Geo.Curve>();
            RH_Geo.Surface surface = null;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref flatPattern)) { return; } ;
            if (!DA.GetData(1, ref initialTarget)) { return; } ;
            if (!DA.GetData(2, ref flatTarget)) { return; } ;

            if (!DA.GetDataList(3, curves)) { return; };
            if (!DA.GetData(4, ref surface)) { return; };

            #region Prepare 

            // ----- Create a Unique List of Ngon and Faces ----- //

            List<RH_Geo.MeshNgon> ngonAndFaces = new List<RH_Geo.MeshNgon>();
            foreach(RH_Geo.MeshNgon ngon in flatPattern.GetNgonAndFacesEnumerable()) { ngonAndFaces.Add(ngon); }

            // ----- Translate Boundary Curves ----- //

            RH_Geo.PolylineCurve[] boundaries = new RH_Geo.PolylineCurve[curves.Count];
            for (int i = 0; i < boundaries.Length; i++) 
            {
                if (curves[i].TryGetPolyline(out RH_Geo.Polyline polyline)) { boundaries[i] = polyline.ToPolylineCurve(); }
                else { throw new ArgumentException("The boundary curves must be polylines.");  }
            }

            #endregion

            #region Create Flat Trimmed Pattern 

            RH_Geo.Mesh flatTrimmedPattern = new RH_Geo.Mesh();

            // ----- Split Edges ----- //

            // Contains the pair of end vertex indices
            // {A} where A = Index of the edge.
            GH.DataTree<(int, int)> i_EdgesInternalSegments = EdgesTrimming(flatPattern, flatTrimmedPattern, boundaries, tolerance, out GH.DataTree<(int, double)> bordersSplitInfos);

            // ----- Face Internal Borders ----- //

            // Contains the consecutive vertex's index composing the face border.
            // {A,B} with A = Index of the face; B = Index of the border relatively to the face. 
            Dictionary<int, List<List<int>>> i_FacesInternalBorders = FacesInternalBorders(i_EdgesInternalSegments, flatPattern, ngonAndFaces, out Dictionary<int, RH_Geo.Polyline> initialFacesBorder);

            // ----- Add Faces ----- //

            foreach (KeyValuePair<int, List<List<int>>> kvp in i_FacesInternalBorders)
            {
                int i_InitialFaceIndex = kvp.Key;
                List<List<int>> i_FaceInternalBorders = kvp.Value;

                for (int i = 0; i < i_FaceInternalBorders.Count; i++)
                {
                    List<int> i_FaceInternalBorder = new List<int>(i_FaceInternalBorders[i]);

                    if (i_FaceInternalBorder[0] == i_FaceInternalBorder[i_FaceInternalBorder.Count - 1]) { i_FaceInternalBorder.RemoveAt(i_FaceInternalBorder.Count - 1); }
                    if (i_FaceInternalBorder.Count < 3) { continue; }

                    AddFace(i_FaceInternalBorder, ref flatTrimmedPattern, globalNormal, tolerance);
                }
            }

            #endregion

            #region Create Trimmed Flat Dual

            // ---------- Create Dual Points ---------- //

            Dictionary<int, RH_Geo.Point3d> facesDualVertices = new Dictionary<int, RH_Geo.Point3d>(i_FacesInternalBorders.Count);
            Dictionary<int, int> facesAssociatedInternalBorder = new Dictionary<int, int>(i_FacesInternalBorders.Count);

            foreach (KeyValuePair<int, List<List<int>>> kvp in i_FacesInternalBorders)
            {
                int i_InitialFaceIndex = kvp.Key;
                List<List<int>> i_FaceInternalBorders = kvp.Value;

                // The face in inside the admissible domain. The dual vertex is the face centroid.
                if (i_FaceInternalBorders.Count == 1 & i_FaceInternalBorders[0][0] == i_FaceInternalBorders[0][i_FaceInternalBorders[0].Count - 1])
                {
                    List<int> i_FaceInternalBorder = i_FaceInternalBorders[0];

                    RH_Geo.Point3d centroid = new RH_Geo.Point3d(0d, 0d, 0d);
                    for (int i = 0; i < i_FaceInternalBorder.Count - 1; i++)
                    {
                        RH_Geo.Point3d position = flatTrimmedPattern.Vertices[i_FaceInternalBorder[i]];
                        centroid += position;
                    }
                    centroid /= i_FaceInternalBorder.Count - 1;

                    facesDualVertices.Add(i_InitialFaceIndex, centroid);
                    facesAssociatedInternalBorder.Add(i_InitialFaceIndex, 0);
                }
                // The face is splitted by the boundaries in one or more borders
                else 
                {
                    RH_Geo.Polyline initialFaceBorder = initialFacesBorder[i_InitialFaceIndex];

                    RH_Geo.Point3d initialFaceCenter = new RH_Geo.Point3d(0d, 0d, 0d);
                    for (int j = 0; j < initialFaceBorder.Count; j++)
                    {
                        initialFaceCenter += initialFaceBorder[j];
                    }
                    initialFaceCenter /= initialFaceBorder.Count;

                    // ---------- Initialisation ---------- //

                    List<int> i_FaceInternalBorder = i_FaceInternalBorders[0];

                    int closestFaceBorderIndex = 0;
                    RH_Geo.Polyline closestFaceBorder = new RH_Geo.Polyline();
                    RH_Geo.Point3d closestFaceCenter = new RH_Geo.Point3d(0d, 0d, 0d);
                    for (int i = 0; i < i_FaceInternalBorder.Count; i++)
                    {
                        RH_Geo.Point3d position = flatTrimmedPattern.Vertices[i_FaceInternalBorder[i]];

                        closestFaceBorder.Add(position);
                        closestFaceCenter += position;
                    }
                    closestFaceBorder.Add(closestFaceBorder[0]);
                    closestFaceCenter /= i_FaceInternalBorder.Count;

                    bool isInside = IsPointInsideCurve(initialFaceCenter, closestFaceBorder.ToPolylineCurve(), tolerance);
                    if (isInside) 
                    { 
                        facesDualVertices.Add(i_InitialFaceIndex, initialFaceCenter);
                        facesAssociatedInternalBorder.Add(i_InitialFaceIndex, 0);
                        continue; 
                    }

                    // ---------- Iteration ---------- //

                    for (int i = 1; i < i_FaceInternalBorders.Count; i++)
                    {
                        List<int> i_FaceOtherInternalBorder = i_FaceInternalBorders[i];

                        RH_Geo.Polyline faceOtherBorder = new RH_Geo.Polyline();
                        RH_Geo.Point3d faceOtherCenter = new RH_Geo.Point3d(0d, 0d, 0d);
                        for (int j = 0; j < i_FaceOtherInternalBorder.Count; j++)
                        {
                            RH_Geo.Point3d position = flatTrimmedPattern.Vertices[i_FaceOtherInternalBorder[j]];

                            faceOtherBorder.Add(position);
                            faceOtherCenter += position;
                        }
                        faceOtherBorder.Add(faceOtherBorder[0]);
                        faceOtherCenter /= i_FaceOtherInternalBorder.Count;

                        isInside = IsPointInsideCurve(initialFaceCenter, faceOtherBorder.ToPolylineCurve(), tolerance);
                        if (isInside) 
                        {
                            facesDualVertices.Add(i_InitialFaceIndex, initialFaceCenter);
                            facesAssociatedInternalBorder.Add(i_InitialFaceIndex, i); 
                            break;
                        }

                        if (faceOtherCenter.DistanceTo(initialFaceCenter) < closestFaceCenter.DistanceTo(initialFaceCenter))
                        {
                            closestFaceCenter = faceOtherCenter; closestFaceBorder = faceOtherBorder; closestFaceBorderIndex = i;
                        }
                    }
                    if (isInside) { continue; }

                    // ---------- Intersect ---------- //

                    RH_Geo.Line line = new RH_Geo.Line(initialFaceCenter, closestFaceCenter - initialFaceCenter);

                    RH_Geo.Point3d closest = closestFaceCenter;

                    RH_Geo.Intersect.CurveIntersections curveIntersections = RH_Geo.Intersect.Intersection.CurveLine(closestFaceBorder.ToPolylineCurve(), line, tolerance, 1e-2 * tolerance);
                    for (int j = 0; j < curveIntersections.Count; j++)
                    {
                        RH_Geo.Point3d projectedPoint = curveIntersections[j].PointA;
                        if (projectedPoint.DistanceTo(initialFaceCenter) < closest.DistanceTo(initialFaceCenter))
                        {
                            closest = projectedPoint;
                        }
                    }

                    facesDualVertices.Add(i_InitialFaceIndex, closest);
                    facesAssociatedInternalBorder.Add(i_InitialFaceIndex, closestFaceBorderIndex);
                }
            }

            // ---------- Heights of Dual Points ---------- //

            List<int> i_InitialFacesIndex = new List<int>(facesDualVertices.Count);
            List<RH_Geo.Point3d> dualVertices = new List<RH_Geo.Point3d>(facesDualVertices.Count);

            foreach (KeyValuePair<int, RH_Geo.Point3d> kvp in facesDualVertices)
            {
                i_InitialFacesIndex.Add(kvp.Key);
                dualVertices.Add(kvp.Value);
            }

            RH_Geo.Brep brep = RH_Geo.Brep.CreateFromSurface(surface);
            RH_Geo.Point3d[] projected = RH_Geo.Intersect.Intersection.ProjectPointsToBreps(new RH_Geo.Brep[1] { brep }, dualVertices, RH_Geo.Vector3d.ZAxis, tolerance);


            /*
                        //Dictionary<int, RH_Geo.Point3d> dualProjectedVertices = new Dictionary<int, RH_Geo.Point3d>(facesDualVertices.Count);
                        //for (int i = 0; i < i_InitialFacesIndex.Length; i++)
                        //{
                        //    dualProjectedVertices.Add(i_InitialFacesIndex[i], projected[i]);
                        //}
            */
            #endregion

            #region Transfer from Flat to Initial

            // ---------- Trimmed Pattern ----------- //

            RH_Geo.Mesh transferedTrimmedPattern = new RH_Geo.Mesh();

            for (int i_V = 0; i_V < flatTrimmedPattern.Vertices.Count; i_V++)
            {
                RH_Geo.Point3d position = flatTrimmedPattern.Vertices[i_V];
                RH_Geo.MeshPoint point = flatTarget.ClosestMeshPoint(position, 1e-2);

                RH_Geo.MeshFace face = initialTarget.Faces[point.FaceIndex];
                if (!face.IsTriangle)
                {
                    this.AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Error, "The barycentric coordinate method is not implemented for none triangular faces.");
                    return;
                }

                RH_Geo.Point3d a = initialTarget.Vertices[face.A];
                RH_Geo.Point3d b = initialTarget.Vertices  [face.B];
                RH_Geo.Point3d c = initialTarget.Vertices[face.C];

                double weightA = point.T[0];
                double weightB = point.T[1];
                double weightC = point.T[2];

                RH_Geo.Point3d transfered = (weightA * a) + (weightB * b) + (weightC * c);

                transferedTrimmedPattern.Vertices.Add(transfered);
            }

            foreach (RH_Geo.MeshFace face in flatTrimmedPattern.Faces) { transferedTrimmedPattern.Faces.AddFace(face); }
            foreach (RH_Geo.MeshNgon ngon in flatTrimmedPattern.Ngons) { transferedTrimmedPattern.Ngons.AddNgon(ngon); }


            // ---------- Trimmed Dual ----------- //

            Dictionary<int, RH_Geo.Point3d> transferedTrimmedDualVertices = new Dictionary<int, RH_Geo.Point3d>(dualVertices.Count);
            for (int i_V = 0; i_V < dualVertices.Count; i_V++)
            {
                RH_Geo.Point3d position = dualVertices[i_V];
                RH_Geo.MeshPoint point = flatTarget.ClosestMeshPoint(position, 0d);

                RH_Geo.MeshFace face = initialTarget.Faces[point.FaceIndex];
                if (!face.IsTriangle)
                {
                    this.AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Error, "The barycentric coordinate method is not implemented for none triangular faces.");
                    return;
                }

                RH_Geo.Point3d a = initialTarget.Vertices[face.A];
                RH_Geo.Point3d b = initialTarget.Vertices[face.B];
                RH_Geo.Point3d c = initialTarget.Vertices[face.C];

                double weightA = point.T[0];
                double weightB = point.T[1];
                double weightC = point.T[2];

                RH_Geo.Point3d transfered = (weightA * a) + (weightB * b) + (weightC * c);

                RH_Geo.Vector3d ab = b - a;
                RH_Geo.Vector3d ac = c - a;
                RH_Geo.Vector3d normal = RH_Geo.Vector3d.CrossProduct(ab, ac);
                normal.Unitize(); 

                transfered += normal * (projected[i_V].Z - dualVertices[i_V].Z);

                transferedTrimmedDualVertices.Add(i_InitialFacesIndex[i_V], transfered); 
            }

            #endregion

            #region Finalise

            // ----- Lines of the Pattern ----- //

            List<RH_Geo.Line> patternLines = new List<RH_Geo.Line>();
            for (int i = 0; i < transferedTrimmedPattern.TopologyEdges.Count; i++)
            {
                RH_Geo.Line edgeLine = transferedTrimmedPattern.TopologyEdges.EdgeLine(i);
                patternLines.Add(edgeLine);
            }
            // ----- Lines of the Dual  ----- //

            He.Mesh<Euc3D.Point> hePrimal = new He.Mesh<Euc3D.Point>();
            Dictionary<int, int> i_OldToNew = new Dictionary<int, int>(flatPattern.Vertices.Count);

            for (int j = 0; j < ngonAndFaces.Count; j++)
            {
                RH_Geo.MeshNgon ngon = ngonAndFaces[j];
                uint[] u_FaceVertices = ngon.BoundaryVertexIndexList();

                List<int> i_FaceVertices = new List<int>(u_FaceVertices.Length);
                for (int i = 0; i < u_FaceVertices.Length; i++)
                {
                    int i_FaceTopoVertex = flatPattern.TopologyVertices.TopologyVertexIndex((int)u_FaceVertices[i]);

                    int i_HeVertex;
                    if (i_OldToNew.ContainsKey(i_FaceTopoVertex)) { i_HeVertex = i_OldToNew[i_FaceTopoVertex]; }
                    else
                    {
                        RH_Geo.Point3d point3d = (RH_Geo.Point3d)flatPattern.TopologyVertices[i_FaceTopoVertex];
                        point3d.CastTo(out Euc3D.Point point);

                        He.Vertex<Euc3D.Point> vertex = hePrimal.AddVertex(point);

                        i_OldToNew.Add(i_FaceTopoVertex, vertex.Index);
                        i_HeVertex = vertex.Index;
                    }
                    i_FaceVertices.Add(i_HeVertex);
                }

                hePrimal.AddFace(i_FaceVertices);
            }

            List<RH_Geo.Line> dualLines = new List<RH_Geo.Line>();
            for (int i = 0; i < hePrimal.FaceCount; i++)
            {
                if (!transferedTrimmedDualVertices.ContainsKey(i)) { continue; }
                RH_Geo.Point3d start = transferedTrimmedDualVertices[i];

                He.Face<Euc3D.Point> face = hePrimal.GetFace(i);
                IReadOnlyList<He.Halfedge<Euc3D.Point>> faceHalfedges = face.FaceHalfedges();
                for (int j = 0; j < faceHalfedges.Count; j++)
                {
                    He.Face<Euc3D.Point> adjacentFace = faceHalfedges[j].PairHalfedge.AdjacentFace;
                    if (adjacentFace == null) { continue; }

                    if (!transferedTrimmedDualVertices.ContainsKey(adjacentFace.Index)) { continue; }

                    RH_Geo.Point3d end = transferedTrimmedDualVertices[adjacentFace.Index];
                    RH_Geo.Line line = new RH_Geo.Line(start, end);
                    dualLines.Add(line);
                }
            }

            // ----- Diagonal Lines ----- //

            GH.DataTree<RH_Geo.Line> diagonalLines = new GH.DataTree<RH_Geo.Line>();
            foreach (KeyValuePair<int, RH_Geo.Point3d> kvp in transferedTrimmedDualVertices)
            {
                int i_InitialFaceIndex = kvp.Key;

                int i_DualFaceBorder = facesAssociatedInternalBorder[i_InitialFaceIndex];
                List<List<int>> i_FaceInternalBorders = i_FacesInternalBorders[i_InitialFaceIndex];

                List<int> i_FaceInternalBorder = i_FaceInternalBorders[i_DualFaceBorder];

                List<RH_Geo.Line> faceLines = new List<RH_Geo.Line>();
                for (int i = 0; i < i_FaceInternalBorder.Count; i++)
                {
                    RH_Geo.Point3d end = transferedTrimmedPattern.Vertices[i_FaceInternalBorder[i]];
                    RH_Geo.Line line = new RH_Geo.Line(kvp.Value, end);
                    faceLines.Add(line);
                }

                diagonalLines.AddRange(faceLines, new GH_Kernel.Data.GH_Path(i_InitialFaceIndex));
            }

            #endregion

            // ----- Set Output ----- //

            DA.SetDataList(0, patternLines);
            DA.SetDataTree(1, diagonalLines);
            DA.SetDataList(2, dualLines);

        }

        #endregion  

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{9E8185EA-E2BD-4A46-94D7-6C40FA5A350C}");

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
        /// <param name="initial"> Mesh being trimmed. </param>
        /// <param name="trimmed"> Mesh resulting from the trim operation. </param>
        /// <param name="boundaries"> Curves defining the boundaries of admissible domain. </param>
        /// <param name="tolerance"> Precision of the computation. </param>
        /// <param name="boundariesSplitInfos"> Information to split the boundaries. (Index of the splitting vertex, Splitting parameter of the border). </param>
        /// <returns> The edges internal segments, usefull for the new trimmed mesh, defined by the indices of the end vertices. </returns>
        private static GH.DataTree<(int, int)> EdgesTrimming(RH_Geo.Mesh initial, RH_Geo.Mesh trimmed, RH_Geo.Curve[] boundaries, double tolerance,
          out GH.DataTree<(int, double)> boundariesSplitInfos)
        {
            GH.DataTree<(int, int)> i_EdgesSegments = new GH.DataTree<ValueTuple<int, int>>();

            RH_Geo.Brep[] breps = RH_Geo.Brep.CreatePlanarBreps(boundaries, tolerance);
            RH_Geo.Brep domain = breps[0];

            // Initialise the necessary datas
            boundariesSplitInfos = new GH.DataTree<ValueTuple<int, double>>();

            Dictionary<int, int> i_OldToNew = new Dictionary<int, int>();

            // Iterate on the edges of the pattern

            int edgeCount = initial.TopologyEdges.Count;
            for (int i = 0; i < edgeCount; i++)
            {
                RH_Geo.Line edgeLine = initial.TopologyEdges.EdgeLine(i);
                RH_Geo.LineCurve edgeCurve = new RH_Geo.LineCurve(edgeLine);

                RH.IndexPair i_EndVertices = initial.TopologyEdges.GetTopologyVertices(i);

                // ----- Compute Intersections ----- //

                List<(int, double, double)> edgeSplitParameters = EdgeSplitParameters(edgeCurve, boundaries, tolerance);

                // ----- Translate Points to Vertex Indices ----- //
                
                List<int> i_EdgeVertices = new List<int>();


                RH_Geo.Point3d projected = domain.ClosestPoint(edgeLine.From);
                bool isStartInside = edgeLine.From.DistanceTo(projected) < tolerance;

                // Manage Start Vertex
                if (isStartInside)
                {
                    int newIndex = ManageExistingVertex(i_EndVertices.I, initial, trimmed, ref i_OldToNew);
                    i_EdgeVertices.Add(newIndex);
                }

                // Manage Intersection Vertices
                foreach ((int, double, double) edgeSplitParameter in edgeSplitParameters)
                {
                    RH_Geo.Point3d position = edgeCurve.PointAt(edgeSplitParameter.Item3);
                    int i_Vertex = trimmed.Vertices.Add(position);

                    i_EdgeVertices.Add(i_Vertex);

                    boundariesSplitInfos.Add((i_Vertex, edgeSplitParameter.Item2), new GH_Kernel.Data.GH_Path(edgeSplitParameter.Item1));
                }

                // Manage End Vertex
                if ((isStartInside & edgeSplitParameters.Count % 2 == 0) | (!isStartInside & edgeSplitParameters.Count % 2 == 1))
                {
                    int newIndex = ManageExistingVertex(i_EndVertices.J, initial, trimmed, ref i_OldToNew);
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
        /// <param name="initial"> Mesh being trimmed. </param>
        /// <param name="trimmed"> Mesh resulting from the trim operation. </param>
        /// <param name="i_OldToNew"> Dictionary defining the conversion of the vertex indices from the initial mesh to the new trimmed mesh. </param>
        /// <returns> The index of the vertex in the new trimmmed mesh. </returns>
        private static int ManageExistingVertex(int i_Vertex, RH_Geo.Mesh initial, RH_Geo.Mesh trimmed, ref Dictionary<int, int> i_OldToNew)
        {
            int i_NewVertex;
            if (i_OldToNew.ContainsKey(i_Vertex))
            {
                i_NewVertex = i_OldToNew[i_Vertex];
            }
            else // Adds the old vertex to the new list
            {
                i_NewVertex = trimmed.Vertices.Add(initial.TopologyVertices[i_Vertex]);
                i_OldToNew.Add(i_Vertex, i_NewVertex);
            }

            return i_NewVertex;
        }

        #endregion

        #region Faces Internal Borders

        /// <summary>
        /// Joins the edge segments of each initial face, and create the border of the trimmed face(s).
        /// </summary>
        /// <param name="i_EdgesInternalSegments"> Edges internal segments usefull for the new trimmed mesh, defined by the indices of the end vertices. </param>
        /// <param name="mesh"> Mesh being trimmed. </param>
        /// <param name="ngonAndFaces"> Ngon and face of the mesh being trimmed. </param>
        /// <param name="initialFacesBorder"> Polyline representing the face border in the initial mesh for the trimmed faces only. </param>
        /// <returns> Indices of the face vertices, forming the internal border of each face. </returns>
        private static Dictionary<int, List<List<int>>> FacesInternalBorders(GH.DataTree<(int, int)> i_EdgesInternalSegments, RH_Geo.Mesh mesh, List<RH_Geo.MeshNgon> ngonAndFaces, 
            out Dictionary<int, RH_Geo.Polyline> initialFacesBorder)
        {
            // Initialise the necessary datas
            Dictionary<int, List<List<int>>> i_FacesInternalBorders = new Dictionary<int, List<List<int>>>();

            //GH.DataTree<int> i_FacesInternalBorders = new GH.DataTree<int>();
            initialFacesBorder = new Dictionary<int, RH_Geo.Polyline>();

            // Iterate on each faces

            for (int i = 0; i < ngonAndFaces.Count; i++)
            {
                RH_Geo.MeshNgon ngon = ngonAndFaces[i];

                //----- Get the ngon edges -----//
                uint[] u_FaceVertices = ngon.BoundaryVertexIndexList();
                int faceVertexCount = u_FaceVertices.Length;

                int[] i_FaceEdges = new int[faceVertexCount];
                RH_Geo.Polyline initialFaceBorder = new RH_Geo.Polyline(faceVertexCount);
                for (int j = 0; j < faceVertexCount; j++)
                {
                    int i_Start = Convert.ToInt32(u_FaceVertices[j]);
                    int i_End = Convert.ToInt32(u_FaceVertices[(j + 1) % faceVertexCount]);

                    int i_TopoStart = mesh.TopologyVertices.TopologyVertexIndex(i_Start);
                    int i_TopoEnd = mesh.TopologyVertices.TopologyVertexIndex(i_End);

                    initialFaceBorder.Add(mesh.TopologyVertices[i_TopoStart]);
                    i_FaceEdges[j] = mesh.TopologyEdges.GetEdgeIndex(i_TopoStart, i_TopoEnd);
                }

                //----- Evaluate whether the edge is flipped. -----//

                bool[] isFaceEdgeFlipped = FaceEdgesOrientation(i_FaceEdges, mesh);

                //----- Combine Face Segments -----//

                List<(int, int)> i_FaceInternalSegments = FaceInternalSegments(i_FaceEdges, isFaceEdgeFlipped, i_EdgesInternalSegments);

                //----- Create Face Borders -----//

                if (i_FaceInternalSegments.Count == 0) { continue; }

                List<List<int>> i_FaceInternalBorders = FaceInternalBorders(i_FaceInternalSegments);

                //----- Stors the Face Borders -----//

                if (1 < i_FaceInternalBorders.Count | (i_FaceInternalBorders.Count == 1 & i_FaceInternalBorders[0][0] != i_FaceInternalBorders[0][i_FaceInternalBorders[0].Count - 1]))
                {
                    initialFacesBorder.Add(i, initialFaceBorder);
                }

                i_FacesInternalBorders.Add(i, i_FaceInternalBorders);
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
        private static List<List<int>> FaceInternalBorders(List<(int, int)> i_FaceInternalSegments)
        {
            List<List<int>> i_FaceInternalBorders = new List<List<int>>();

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
                    i_FaceInternalBorders.Add(i_FaceInternalBorder);

                    i_FaceInternalBorder = new List<int> { i_FaceSegment.Item1, i_FaceSegment.Item2 };
                }
            }
            i_FaceInternalBorders.Add(i_FaceInternalBorder);

            // Finalise : Joins the last and the first borders if needed
            if (i_FaceInternalBorders.Count > 1)
            {
                GH_Kernel.Data.GH_Path firstPath = new GH_Kernel.Data.GH_Path(0);
                GH_Kernel.Data.GH_Path lastPath = new GH_Kernel.Data.GH_Path(i_FaceInternalBorders.Count - 1);

                List<int> first = i_FaceInternalBorders[0];
                List<int> last = i_FaceInternalBorders[i_FaceInternalBorders.Count - 1];

                if (first[0] == last[last.Count - 1])
                {
                    first.RemoveAt(0);
                    last.AddRange(first);
                    i_FaceInternalBorders.RemoveAt(0);
                }
            }

            return i_FaceInternalBorders;
        }

        #endregion

        #region Add Face

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

        #endregion
    }

}
