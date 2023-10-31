using System;
using System.Collections.Generic;

using Euc3D = BRIDGES.Geometry.Euclidean3D;
using He = BRIDGES.DataStructures.PolyhedralMeshes.HalfedgeMesh;
using GP = BRIDGES.Solvers.GuidedProjection;
using LinAlg_Vect = BRIDGES.LinearAlgebra.Vectors;

using RH_Geo = Rhino.Geometry;

using GH_Kernel = Grasshopper.Kernel;

using BRIDGES.McNeel.Rhino.Extensions.Geometry.Euclidean3D;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;
using Gh_Param_Euc3D = BRIDGES.McNeel.Grasshopper.Parameters.Geometry.Euclidean3D;
using Gh_Types_Euc3D = BRIDGES.McNeel.Grasshopper.Types.Geometry.Euclidean3D;

using Types_Var = Llama.Types.Variables;
using Types_En = Llama.Types.Energies;
using Params_Var = Llama.Parameters.Variables;
using Params_En = Llama.Parameters.Energies;


namespace Llama.Energies
{
    #region Component Hexagonal Mesh

    /// <summary>
    /// A grasshopper component giving the fan orientation in a mesh.
    /// </summary>
    internal class Comp_HexagonalMesh : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_FanOrientation"/> class.
        /// </summary>
        public Comp_HexagonalMesh()
          : base("Hexagonal Mesh", "Hexa",
              "Creates an hexagonal halfedge mesh from faces as polylines.",
              Settings.CategoryName, "Nexorade")
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polylines", "P", "Polylines representing the faces of the hewagonal mesh", GH_Kernel.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Gh_Param_Euc3D.Param_HeMesh(), "Mesh", "M", "Halfedge mesh.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/
            List<RH_Geo.Curve> curves = new List<RH_Geo.Curve>();

            /******************** Get Inputs ********************/

            if (!DA.GetDataList(0, curves)) { return; };

            /******************** Core ********************/

            He.Mesh<Euc3D.Point> mesh = new He.Mesh<Euc3D.Point>();

            foreach (RH_Geo.Curve curve in curves)
            {
                curve.TryGetPolyline(out RH_Geo.Polyline polyline);

                int vertexCount = polyline.IsClosed ? polyline.Count - 1 : polyline.Count;

                List<He.Vertex<Euc3D.Point>> faceVertices = new List<He.Vertex<Euc3D.Point>>(vertexCount);
                for (int i = 0; i < vertexCount; i++)
                {
                    polyline[i].CastTo(out Euc3D.Point point);

                    He.Vertex<Euc3D.Point> faceVertex = null;
                    foreach (He.Vertex<Euc3D.Point> vertex in mesh.GetVertices())
                    {
                        if (vertex.Position.DistanceTo(point) < 1e-4) { faceVertex = vertex; break; }
                    }

                    if (faceVertex is null) { faceVertex = mesh.AddVertex(point); }

                    faceVertices.Add(faceVertex);
                }

                mesh.AddFace(faceVertices);

            }

            /******************** Set Output ********************/

            DA.SetData(0, mesh);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{AA7DDB73-BC5D-4BB8-A7BA-A42AECB2A8B2}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure.primary);

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

    #endregion

    #region Component Fan Orientation

    /// <summary>
    /// A grasshopper component giving the fan orientation in a mesh.
    /// </summary>
    public class Comp_FanOrientation : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_FanOrientation"/> class.
        /// </summary>
        public Comp_FanOrientation()
          : base("Fan Orientation", "Fan Ori.",
              "Gives the orientation of the fans in a mesh.",
              Settings.CategoryName, "Nexorade")
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Gh_Param_Euc3D.Param_HeMesh(), "Mesh", "M", "Halfedge mesh.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Fan Orientation", "O", "Fan Orientation : 1 = Left; -1 = Right; 0 = Not Assigned. All the fans should be assigned.", GH_Kernel.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/

            Gh_Types_Euc3D.Gh_HeMesh gh_Mesh = null;

            /******************** Get Inputs ********************/

            if (!DA.GetData(0, ref gh_Mesh)) { return; };
            He.Mesh<Euc3D.Point> heMesh = gh_Mesh.Value;

            /******************** Core ********************/

            int vertexCount = heMesh.VertexCount;

            // Fan Orientation : 1 = Left; -1 = Right; 0 = Not Assigned
            int[] fanOrientation = new int[vertexCount];

            foreach (He.Vertex<Euc3D.Point> v in heMesh.GetVertices())
            {
                if (v.OutgoingHalfedge is null) { heMesh.RemoveVertex(v); }
            }

            He.Vertex<Euc3D.Point> v_First = heMesh.GetVertex(0);
            He.Halfedge<Euc3D.Point> he_First = v_First.OutgoingHalfedge;
            fanOrientation[he_First.StartVertex.Index] = 1;
            fanOrientation[he_First.EndVertex.Index] = -1;

            SpreadOrientation(he_First);

            // Halfedge whose:
            // - Start Vertex is being used for propagation (turning around it)
            // - End Vertex is where the propagation comes from (i.e. its fan is already assigned)
            void SpreadOrientation(He.Halfedge<Euc3D.Point> halfedge)
            {
                int orientation = fanOrientation[halfedge.EndVertex.Index];

                He.Halfedge<Euc3D.Point> he = halfedge.PrevHalfedge.PairHalfedge;

                while (he != halfedge)
                {
                    if (fanOrientation[he.EndVertex.Index] == 0)
                    {
                        fanOrientation[he.EndVertex.Index] = orientation;
                        SpreadOrientation(he.PairHalfedge);
                    }

                    he = he.PrevHalfedge.PairHalfedge;
                }
            }

            /******************** Set Output ********************/

            DA.SetDataList(0, fanOrientation);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{FE32D8D4-F03F-4CAC-B651-3D906B64DD36}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => GH_Kernel.GH_Exposure.primary;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

    #endregion

    #region Component Nexorade

    /// <summary>
    /// A grasshopper component building the energies for the nexorade.
    /// </summary>
    public class Comp_Nexorade : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_Nexorade"/> class.
        /// </summary>
        public Comp_Nexorade()
          : base("Nexorade", "Nexorade",
              "Creates the Nexorade Energies for the GPA Solver.",
              Settings.CategoryName, "Nexorade")
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Gh_Param_Euc3D.Param_HeMesh(), "Mesh", "M", "Halfedge mesh.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Fan Orientation", "O", "Fan Orientation : 1 = Left; -1 = Right; 0 = Not Assigned. All the fans should be assigned.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddNumberParameter("Eccentricity", "E", "Eccentricity value for the fans.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Engagement Length", "L", "Engagement Length value for the fans.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params_Var.Param_VariableSet(), "Set", "S", "Variable set containing the translation vectors of each edges", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params_En.Param_Energy(), "Energy", "E", "Energies to control the fans.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialisation ----- //

            Gh_Types_Euc3D.Gh_HeMesh gh_HeMesh = null;
            List<int> fanOrientation = new List<int>();
            double eccentricity = 0d;
            double engagementLength = 0d;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref gh_HeMesh)) { return; };
            He.Mesh<Euc3D.Point> heMesh = gh_HeMesh.Value;

            if (!DA.GetDataList(1, fanOrientation)) { return; };
            if (!DA.GetData(2, ref eccentricity)) { return; };
            if (!DA.GetData(3, ref engagementLength)) { return; };

            // ----- Prepare ----- //

            // Translation vector for each edge of the halfedge mesh
            GP.Variable[] translations = new GP.Variable[heMesh.EdgeCount];
            for (int i = 0; i < translations.Length; i++) { translations[i] = new GP.Variable(0d, 0d, 0d); }

            // ----- Energies ----- //

            List<GP.Energy> energies = new List<GP.Energy>(heMesh.EdgeCount + 2 * heMesh.HalfedgeCount);

            #region Reduce Space

            for (int i = 0; i < heMesh.EdgeCount; i++)
            {
                He.Halfedge<Euc3D.Point> he = heMesh.GetHalfedge(2 * i);
                Euc3D.Vector vector = he.EndVertex.Position - he.StartVertex.Position;
                vector.Unitise();

                double[] components = new double[3] { vector.X, vector.Y, vector.Z };
                VectorOrthogonality energyType = new VectorOrthogonality(components);

                GP.Energy energy = new GP.Energy(energyType, new GP.Variable[1] { translations[i] }, 10d);
                energies.Add(energy);
            }

            #endregion

            #region Regularisation

            for (int i = 0; i < heMesh.EdgeCount; i++)
            {
                He.Halfedge<Euc3D.Point> he = heMesh.GetHalfedge(2 * i);

                Euc3D.Vector normal = new Euc3D.Vector(0d, 0d, 0d);
                if (he.AdjacentFace != null) { normal += FaceNormal(he.AdjacentFace); }
                if (he.PairHalfedge.AdjacentFace != null) { normal += FaceNormal(he.PairHalfedge.AdjacentFace); }
                normal.Unitise();

                double[] components = new double[3] { normal.X, normal.Y, normal.Z };
                VectorOrthogonality energyType = new VectorOrthogonality(components);

                GP.Energy energy = new GP.Energy(energyType, new GP.Variable[1] { translations[i] }, 1d);
                energies.Add(energy);

                // Index of the face to get the normal
                Euc3D.Vector FaceNormal(He.Face<Euc3D.Point> face)
                {
                    IReadOnlyList<He.Vertex<Euc3D.Point>> faceVertices = face.FaceVertices();

                    Euc3D.Point faceCenter = new Euc3D.Point(0d, 0d, 0d);
                    for (int j = 0; j < faceVertices.Count; j++) { faceCenter += faceVertices[j].Position; }
                    faceCenter /= faceVertices.Count;

                    Euc3D.Vector faceNormal = new Euc3D.Vector(0d, 0d, 0d);
                    for (int j = 0; j < faceVertices.Count; j++)
                    {
                        int j1 = (j + 1) % faceVertices.Count;

                        Euc3D.Vector vector = Euc3D.Vector.CrossProduct(faceVertices[j].Position - faceCenter, faceVertices[j1].Position - faceCenter);
                        vector.Unitise();
                        faceNormal += vector;
                    }

                    faceNormal /= faceVertices.Count;
                    return faceNormal;
                }

            }

            #endregion

            #region Eccentricity Energy

            foreach (He.Halfedge<Euc3D.Point> he in heMesh.GetHalfedges())
            {
                He.Vertex<Euc3D.Point> start = he.StartVertex;

                He.Halfedge<Euc3D.Point> leftHe = he.PairHalfedge.NextHalfedge;
                He.Halfedge<Euc3D.Point> rightHe = he.PrevHalfedge.PairHalfedge;

                He.Halfedge<Euc3D.Point> otherHe = leftHe;

                double localEccentricity;
                if (fanOrientation[start.Index] == 0) { throw new Exception("The fan Orientation must be encoded by either 1 (Left) or -1 (Right)."); }
                else { localEccentricity = eccentricity * fanOrientation[start.Index]; }

                if ((he.IsBoundary() | he.PairHalfedge.IsBoundary()) & (otherHe.IsBoundary() | otherHe.PairHalfedge.IsBoundary()))
                {
                    /* There is no need to create an energy in this configuration */
                }
                else
                {
                    Euc3D.Vector vector = he.EndVertex.Position - he.StartVertex.Position;
                    Euc3D.Vector otherVector = otherHe.EndVertex.Position - otherHe.StartVertex.Position;

                    int edgeIndex = he.Index % 2 == 0 ? he.Index / 2 : (he.Index - 1) / 2;
                    int otherEdgeIndex = otherHe.Index % 2 == 0 ? otherHe.Index / 2 : (otherHe.Index - 1) / 2;

                    Eccentricity energyType = new Eccentricity(vector, otherVector, localEccentricity);

                    GP.Energy energy = new GP.Energy(energyType, new GP.Variable[2] { translations[edgeIndex], translations[otherEdgeIndex] }, 50d);
                    energies.Add(energy);
                }
            }

            #endregion

            #region Engagement Length Energy

            foreach (He.Halfedge<Euc3D.Point> he in heMesh.GetHalfedges())
            {
                He.Vertex<Euc3D.Point> start = he.StartVertex;

                He.Halfedge<Euc3D.Point> leftHe = he.PairHalfedge.NextHalfedge;
                He.Halfedge<Euc3D.Point> rightHe = he.PrevHalfedge.PairHalfedge;

                if (start.Valence() < 3 | he.IsBoundary() | he.PairHalfedge.IsBoundary())
                {
                    /* The engagement Length can't be computed with only two edges.
                     * The energy should not be expressed on the boundaries */
                }
                else
                {
                    Euc3D.Vector leftVector = leftHe.EndVertex.Position - leftHe.StartVertex.Position;
                    Euc3D.Vector middleVector = he.EndVertex.Position - he.StartVertex.Position;
                    Euc3D.Vector rightVector = rightHe.EndVertex.Position - rightHe.StartVertex.Position;

                    if (fanOrientation[start.Index] == -1 /* Right Fan ? */) { leftVector = -leftVector; middleVector = -middleVector; rightVector = -rightVector; }

                    int l = leftHe.Index % 2 == 0 ? leftHe.Index / 2 : (leftHe.Index - 1) / 2;
                    int m = he.Index % 2 == 0 ? he.Index / 2 : (he.Index - 1) / 2;
                    int r = rightHe.Index % 2 == 0 ? rightHe.Index / 2 : (rightHe.Index - 1) / 2;

                    EngagementLength energyType = new EngagementLength(leftVector, middleVector, rightVector, engagementLength);

                    GP.Energy energy = new GP.Energy(energyType, new GP.Variable[3] { translations[l], translations[m], translations[r] }, 50d);
                    energies.Add(energy);
                }

            }

            #endregion

            // ----- Set Output ----- //

            Types_Var.Gh_VariableSet set = new Types_Var.Gh_VariableSet(translations, "Translations");

            DA.SetData(0, set);

            DA.SetDataList(1, energies);


            //DA.SetDataList(0, fanOrientation);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{0DF18FA0-FF3C-4710-AA7B-5A92A1C5A40F}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => GH_Kernel.GH_Exposure.secondary;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

    #endregion

    #region Component Rebuild

    /// <summary>
    /// A grasshopper component rebuilding the nexorade.
    /// </summary>
    public class Comp_Rebuild : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_Rebuild"/> class.
        /// </summary>
        public Comp_Rebuild()
          : base("Rebuild after Solving", "Rebuild",
              "Rebuilds the nexorade from the translation vectors.",
              Settings.CategoryName, "Nexorade")
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Gh_Param_Euc3D.Param_HeMesh(), "Mesh", "M", "Halfedge mesh.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Fan Orientation", "O", "Fan Orientation : 1 = Left; -1 = Right; 0 = Not Assigned. All the fans should be assigned.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddNumberParameter("Variables Components", "V", "Components of the variables representing the translation vectors of each edges", GH_Kernel.GH_ParamAccess.tree);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Lines", "L", "Lines of the Nexorade", GH_Kernel.GH_ParamAccess.list);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/

            // ----- Initialisation ----- //

            Gh_Types_Euc3D.Gh_HeMesh gh_Mesh = null;
            List<int> fanOrientation = new List<int>();
            GH_Kernel.Data.GH_Structure<GH_Kernel.Types.GH_Number> gh_Structure = null;


            /******************** Get Inputs ********************/

            if (!DA.GetData(0, ref gh_Mesh)) { return; };
            He.Mesh<Euc3D.Point> heMesh = gh_Mesh.Value;

            if (!DA.GetDataList(1, fanOrientation)) { return; };
            if (!DA.GetDataTree(2, out gh_Structure)) { return; };

            /******************** Core ********************/

            Euc3D.Vector[] translations = new Euc3D.Vector[gh_Structure.PathCount];
            for (int i = 0; i < gh_Structure.PathCount; i++)
            {
                List<GH_Kernel.Types.GH_Number> branch = gh_Structure[i];

                translations[i] = new Euc3D.Vector(branch[0].Value, branch[1].Value, branch[2].Value);
            }

            List<RH_Geo.Line> lines = new List<RH_Geo.Line>();
            for (int index = 0; index < heMesh.EdgeCount; index++)
            {
                He.Halfedge<Euc3D.Point> halfedge = heMesh.GetHalfedge(2 * index);

                if (halfedge.IsBoundary() | halfedge.PairHalfedge.IsBoundary()) { continue; }

                Euc3D.Point start = ComputeNewStart(2 * index);
                Euc3D.Point end = ComputeNewStart(2 * index + 1);

                start.CastTo(out RH_Geo.Point3d rh_Start);
                end.CastTo(out RH_Geo.Point3d rh_End);

                RH_Geo.Line line = new RH_Geo.Line(rh_Start, rh_End);
                lines.Add(line);
            }

            // Halfedge Index 
            Euc3D.Point ComputeNewStart(int index)
            {
                He.Halfedge<Euc3D.Point> he = heMesh.GetHalfedge(index);
                Euc3D.Vector vector = he.EndVertex.Position - he.StartVertex.Position;
                vector.Unitise();

                He.Halfedge<Euc3D.Point> otherHe;
                if (fanOrientation[he.StartVertex.Index] == 1 /* Left Fan ? */) { otherHe = he.PrevHalfedge.PairHalfedge; }
                else if (fanOrientation[he.StartVertex.Index] == -1 /* Right Fan ? */) { otherHe = he.PairHalfedge.NextHalfedge; }
                else { throw new Exception("The fan orientation must be encoded by either 1 (Left) or -1 (Right)."); }

                Euc3D.Vector otherVector = otherHe.EndVertex.Position - otherHe.StartVertex.Position;

                int edgeIndex = he.Index % 2 == 0 ? he.Index / 2 : (he.Index - 1) / 2;
                int otherEdgeIndex = otherHe.Index % 2 == 0 ? otherHe.Index / 2 : (otherHe.Index - 1) / 2;

                Euc3D.Point heStart = he.StartVertex.Position;

                Euc3D.Vector translation = translations[edgeIndex];
                Euc3D.Vector otherTranslation = translations[otherEdgeIndex];

                // Works great for orthogonal direction but not otherwise ? 
                // double factor = Euc3D.Vector.DotProduct(otherTranslation - translation, vector);

                // From non-orthogonal projection

                if ((he.IsBoundary() | he.PairHalfedge.IsBoundary()) & (otherHe.IsBoundary() | otherHe.PairHalfedge.IsBoundary()))
                {
                    return heStart;
                }
                else
                {
                    Euc3D.Vector deltaT = -(translation - otherTranslation);

                    double factor = (deltaT * vector * (otherVector * otherVector) - deltaT * otherVector * (vector * otherVector))
                        / (vector * vector * (otherVector * otherVector) - vector * otherVector * (vector * otherVector));

                    return heStart + translation + factor * vector;
                }
            }

            /******************** Set Output ********************/

            DA.SetDataList(0, lines);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{DB7A1839-48D3-498D-8B23-F65CB10AAB00}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => GH_Kernel.GH_Exposure.tertiary;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

    #endregion


    #region Energies & Constraints

    /// <summary>
    /// Energy enforcing two initially co-planar three-dimensional euclidean ray to have a fixed eccentricity. The list of variables of this energy consists in:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>T<sub>i</sub></term>
    ///         <description> Variable representing the translation vector of the first ray.</description>
    ///     </item>
    ///     <item> 
    ///         <term>T<sub>j</sub></term>
    ///         <description> Variable representing the translation vector of the second ray.</description>
    ///     </item>
    /// </list>
    /// </summary>
    public class Eccentricity : GP.Abstracts.EnergyType
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Eccentricity"/> class.
        /// </summary>
        /// <param name="axis">  Axis of the first ray.</param>
        /// <param name="otherAxis"> Axis of the second ray. </param>
        /// <param name="eccentricity"> Target eccentricity between the three-dimensional euclidean ray. </param>
        /// <remarks> The ray axes should be "outgoing" ... </remarks>
        public Eccentricity(Euc3D.Vector axis, Euc3D.Vector otherAxis, double eccentricity)
        {
            // ----- Define Ki ----- //

            int[] rowIndices = new int[6] { 0, 1, 2, 3, 4, 5 };

            Euc3D.Vector crossProduct = Euc3D.Vector.CrossProduct(axis, otherAxis); crossProduct.Unitise();
            double[] values = new double[6]
            {
                - crossProduct.X, - crossProduct.Y, - crossProduct.Z,
                crossProduct.X, crossProduct.Y, crossProduct.Z,
            };

            LocalKi = new LinAlg_Vect.SparseVector(6, rowIndices, values);

            // ----- Define Si ----- //

            Si = eccentricity;
        }

        #endregion
    }

    /// <summary>
    /// Energy enforcing two initially co-planar three-dimensional euclidean ray to have a fixed length along a third middle ray (also contained in the inital plane). The list of variables of this energy consists in:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>T<sub>l</sub></term>
    ///         <description> Variable representing the translation vector of the left ray.</description>
    ///     </item>
    ///     <item> 
    ///         <term>T<sub>m</sub></term>
    ///         <description> Variable representing the translation vector of the middle ray.</description>
    ///     </item>
    ///     <item> 
    ///         <term>T<sub>r</sub></term>
    ///         <description> Variable representing the translation vector of the right ray.</description>
    ///     </item>
    /// </list>
    /// </summary>
    public class EngagementLength : GP.Abstracts.EnergyType
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Eccentricity"/> class.
        /// </summary>
        /// <param name="leftAxis"> Axis of the left ray. </param>
        /// <param name="middleAxis"> Axis of the middle ray. </param>
        /// <param name="rightAxis">  Axis of the right ray.</param>
        /// <param name="length"> Target engagement length.. </param>
        /// <remarks> For left fans, the rays should represent ougoing directions. Otherwise, the ray axes should be fliped for right fans. </remarks>
        public EngagementLength(Euc3D.Vector leftAxis, Euc3D.Vector middleAxis, Euc3D.Vector rightAxis, double length)
        {
            // ----- Define Ki ----- //

            int[] rowIndices = new int[9] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

            leftAxis.Unitise(); middleAxis.Unitise(); rightAxis.Unitise();
            double[] values = new double[9];

            double sl_Left = leftAxis * leftAxis;
            double sl_Middle = middleAxis * middleAxis, l_Middle = Math.Sqrt(sl_Middle);
            double sl_Right = rightAxis * rightAxis;

            double dot_MiddleLeft = middleAxis * leftAxis;
            double dot_MiddleRight = middleAxis * rightAxis;

            double d_MiddleLeft = sl_Middle * sl_Left - dot_MiddleLeft * dot_MiddleLeft * (dot_MiddleLeft * dot_MiddleLeft);
            //d_MiddleLeft /= l_Middle;
            double d_MiddleRight = sl_Middle * sl_Right - dot_MiddleRight * dot_MiddleRight * (dot_MiddleRight * dot_MiddleRight);
            //d_MiddleRight /= l_Middle;

            double f_Left = sl_Left / d_MiddleLeft;
            double f_MiddleLeft = dot_MiddleLeft / d_MiddleLeft;
            double f_MiddleRight = dot_MiddleRight / d_MiddleRight;
            double f_Right = sl_Right / d_MiddleRight;

            // For X coordinates
            double tmp_MiddleLeft = f_Left * middleAxis.X - f_MiddleLeft * leftAxis.X;
            double tmp_MiddleRight = f_Right * middleAxis.X - f_MiddleRight * rightAxis.X;
            values[0] = tmp_MiddleLeft;
            values[3] = -tmp_MiddleLeft + tmp_MiddleRight;
            values[6] = -tmp_MiddleRight;

            // For Y coordinates
            tmp_MiddleLeft = f_Left * middleAxis.Y - f_MiddleLeft * leftAxis.Y;
            tmp_MiddleRight = f_Right * middleAxis.Y - f_MiddleRight * rightAxis.Y;
            values[1] = tmp_MiddleLeft;
            values[4] = -tmp_MiddleLeft + tmp_MiddleRight;
            values[7] = -tmp_MiddleRight;

            // For Z coordinates
            tmp_MiddleLeft = f_Left * middleAxis.Z - f_MiddleLeft * leftAxis.Z;
            tmp_MiddleRight = f_Right * middleAxis.Z - f_MiddleRight * rightAxis.Z;
            values[2] = tmp_MiddleLeft;
            values[5] = -tmp_MiddleLeft + tmp_MiddleRight;
            values[8] = -tmp_MiddleRight;

            LocalKi = new LinAlg_Vect.SparseVector(9, rowIndices, values);

            // ----- Define Si ----- //

            Si = length;
        }

        #endregion
    }

    /// <summary>
    /// Energy enforcing a vector variable to be orthogonal to a fixed direction. The list of variables of this energy consists in:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>V</term>
    ///         <description> Variable representing the vector.</description>
    ///     </item>
    /// </list>
    /// </summary>
    public class VectorOrthogonality : GP.Abstracts.EnergyType
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="VectorOrthogonality"/> class.
        /// </summary>
        /// <param name="direction"> Components of the direction. </param>
        public VectorOrthogonality(double[] direction)
        {
            // ----- Define Ki ----- //

            int[] rowIndices = new int[direction.Length];
            double[] values = new double[direction.Length];
            for (int i = 0; i < direction.Length; i++)
            {
                rowIndices[i] = i; values[i] = direction[i];
            }

            LocalKi = new LinAlg_Vect.SparseVector(direction.Length, rowIndices, values);

            // ----- Define Si ----- //

            Si = 0d;
        }

        #endregion
    }

    #endregion
}
