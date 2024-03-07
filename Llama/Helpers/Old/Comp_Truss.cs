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
    internal class Comp_Truss : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_Truss"/> class.
        /// </summary>
        public Comp_Truss()
          : base("Truss", "Truss",
              "Truss the primal and dual pattern together.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Helpers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Primal Vertices", "Vp", "Vertices of the flat trimmed pattern.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddIntegerParameter("Primal Faces", "Fp", "Faces of the flat trimmed pattern.", GH_Kernel.GH_ParamAccess.tree);

            pManager.AddPointParameter("Dual Vertices", "Vd", "Vertices of the dual.", GH_Kernel.GH_ParamAccess.tree);
            pManager.AddPathParameter("Dual Faces", "Fd", "Faces vertices of the dual.", GH_Kernel.GH_ParamAccess.tree);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Transfered Primal", "Pt", "Primal Pattern transfered onto the target's initial configuration.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddLineParameter("Diagonals", "D", "Diagonals connecting the Primal Pattern with the Dual Pattern", GH_Kernel.GH_ParamAccess.tree);
            pManager.AddMeshParameter("Transfered Dual", "Dt", "Dual Pattern transfered onto the target's initial configuration with offset.", GH_Kernel.GH_ParamAccess.item);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            double tolerance = 1e-4;


            // ----- Initialise ----- //

            List<RH_Geo.Point3d> primalVertices = new List<RH_Geo.Point3d>();
            GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Integer> struct_PrimalFaces = new GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Integer>();

            GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Point> struct_DualVertices = new GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Point>();
            GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_StructurePath> struct_DualFaces = new GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_StructurePath>();


            // ----- Get Inputs ----- //

            if (!DA.GetDataList(0, primalVertices)) { return; } ;
            if (!DA.GetDataTree(1, out struct_PrimalFaces)) { return; } ;

            if (!DA.GetDataTree(2, out struct_DualVertices)) { return; } ;
            if (!DA.GetDataTree(3, out struct_DualFaces)) { return; } ;

            #region GH_Structure to DataTree

            GH.DataTree<int> i_PrimalFacesVertices = new GH.DataTree<int>();
            for (int i = 0; i < struct_PrimalFaces.PathCount; i++)
            {
                GH_Kernel.Data.GH_Path path = struct_PrimalFaces.Paths[i];
                List<GH_Kernel.Types.GH_Integer> branch = struct_PrimalFaces[path];
                for (int j = 0; j < branch.Count; j++)
                {
                    i_PrimalFacesVertices.Add(branch[j].Value, path);
                }
            }

            GH.DataTree<RH_Geo.Point3d> dualVertices = new GH.DataTree<RH_Geo.Point3d>();
            for (int i = 0; i < struct_DualVertices.PathCount; i++)
            {
                GH_Kernel.Data.GH_Path path = struct_DualVertices.Paths[i];
                List<GH_Kernel.Types.GH_Point> branch = struct_DualVertices[path];
                for (int j = 0; j < branch.Count; j++)
                {
                    dualVertices.Add(branch[j].Value, path);
                }
            }

            GH.DataTree<GH_Kernel.Data.GH_Path> p_DualFacesVertices = new GH.DataTree<GH_Kernel.Data.GH_Path>();
            for (int i = 0; i < struct_DualFaces.PathCount; i++)
            {
                GH_Kernel.Data.GH_Path path = struct_DualFaces.Paths[i];
                List<GH_Kernel.Types.GH_StructurePath> branch = struct_DualFaces[path];
                for (int j = 0; j < branch.Count; j++)
                {
                    p_DualFacesVertices.Add(branch[j].Value, path);
                }
            }

            #endregion

            // ----- Core ----- //

            #region Create Primal Meshes

            // ----- Flat Primal Mesh ----- //

            RH_Geo.Mesh primal = new RH_Geo.Mesh();
            for (int i = 0; i < primalVertices.Count; i++)
            {
                primal.Vertices.Add(primalVertices[i]);
            }
            for (int i = 0; i < i_PrimalFacesVertices.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = i_PrimalFacesVertices.Paths[i];

                List<int> i_PrimalFaceVertices = i_PrimalFacesVertices.Branch(path);

                AddFace(i_PrimalFaceVertices, ref primal, tolerance);
            }

            #endregion

            #region Create Transfered Dual Mesh

            // ----- Flat Dual Mesh ----- //

            RH_Geo.Mesh dual = new RH_Geo.Mesh();

            Dictionary<GH_Kernel.Data.GH_Path, int> i_InitialFaceToDualVertex = new Dictionary<GH_Kernel.Data.GH_Path, int>();

            GH.DataTree<int> tree = new GH.DataTree<int>();
            for (int i = 0; i < dualVertices.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = dualVertices.Paths[i];

                RH_Geo.Point3d position = dualVertices[path, 0];

                int i_Vertex = dual.Vertices.Add(position);
                i_InitialFaceToDualVertex.Add(path, i_Vertex);
                tree.Add(i_Vertex, path);
            }
            for (int i = 0; i < p_DualFacesVertices.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = p_DualFacesVertices.Paths[i];

                List<GH_Kernel.Data.GH_Path> p_DualFace = p_DualFacesVertices.Branch(path);

                List<int> i_DualFaceVertices = new List<int>();
                for (int j = 0; j < p_DualFace.Count; j++)
                {
                    //i_DualFaceVertices.Add(i_InitialFaceToDualVertex[p_DualFace[j]]);
                    i_DualFaceVertices.Add(tree[p_DualFace[j], 0]);
                }

                AddFace(i_DualFaceVertices, ref dual, tolerance);
            }

            #endregion

            #region Create Diagonals


            GH.DataTree<RH_Geo.Line> diagonals = new GH.DataTree<RH_Geo.Line>();
            for (int i = 0; i < i_PrimalFacesVertices.Paths.Count; i++)
            {
                GH_Kernel.Data.GH_Path path = i_PrimalFacesVertices.Paths[i];

                RH_Geo.Point3d dualPoint = dualVertices[path, 0];

                List<int> i_PrimalFaceVertices = i_PrimalFacesVertices.Branch(path);
                for (int j = 0; j < i_PrimalFaceVertices.Count; j++)
                {
                    RH_Geo.Point3d primalPoint = primalVertices[i_PrimalFaceVertices[j]];

                    RH_Geo.Line line = new RH_Geo.Line(dualPoint, primalPoint);

                    diagonals.Add(line, path);
                }
            }

            #endregion

            // ----- Set Output ----- //

            DA.SetData(0, primal);
            DA.SetDataTree(1, diagonals);
            DA.SetData(2, dual);

        }

        #endregion  

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{B81158C1-8C9F-430A-A4AF-6326A1EF722E}");

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


        #region Add Face

        /// <summary>
        /// Add a face to the mesh from the vertex indices.
        /// </summary>
        /// <param name="mesh"> Mesh in which the face should be added. </param>
        /// <param name="i_FaceVertices"> Vertex indices of the face to add. </param>
        /// <param name="tolerance"> Brep representing the admissible domain for the pattern. </param>
        void AddFace(List<int> i_FaceVertices, ref RH_Geo.Mesh mesh, double tolerance)
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
                            RH_Geo.Intersect.CurveIntersections curveIntersections = RH_Geo.Intersect.Intersection.CurveCurve(border, lineCurve, tolerance,  tolerance);
                            foreach (RH_Geo.Intersect.IntersectionEvent intersection in curveIntersections)
                            {
                                if (intersection.IsOverlap)
                                {
                                    if ((tolerance < intersection.PointA.DistanceTo(first) & tolerance < intersection.PointA.DistanceTo(last))
                                        | (tolerance < intersection.PointB.DistanceTo(first) & tolerance < intersection.PointB.DistanceTo(last)))
                                    {
                                        isNotSplitted = false; break;
                                    }
                                }
                                else if (intersection.IsPoint && (tolerance < intersection.PointA.DistanceToSquared(first) & tolerance < intersection.PointA.DistanceToSquared(last)))
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

            bool isWellOriented = true;

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
