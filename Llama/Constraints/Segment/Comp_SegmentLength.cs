using System;

using GP = BRIDGES.Solvers.GuidedProjection;
using LinAlg_Mat = BRIDGES.LinearAlgebra.Matrices;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Types_Var = Llama.Types.Variables;
using Types_Cstr = Llama.Types.Constraints;
using Params_Var = Llama.Parameters.Variables;
using Params_Cstr = Llama.Parameters.Constraints;


namespace Llama.Constraints.Segment
{
    /// <summary>
    /// A grasshopper component creating a <see cref="SegmentLength"/>-based <see cref="GP.Constraint"/>..
    /// </summary>
    public class Comp_SegmentLength : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_SegmentLength"/> class.
        /// </summary>
        public Comp_SegmentLength()
          : base("Segment Length", "Length",
              "Create a Segment Length constraint for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Constraints])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params_Var.Param_Variable(), "Start Variable", "S", "Variable representing the start of the segment", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params_Var.Param_Variable(), "End Variable", "E", "Variable representing the end of the segment", GH_Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Taget Length", "L", "Target length of the segment.", GH_Kernel.GH_ParamAccess.item);

            pManager.AddNumberParameter("Weight", "W", "Weight of the constraint.", GH_Kernel.GH_ParamAccess.item);

            pManager[3].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params_Cstr.Param_Constraint(), "Constraint", "C", "Segment Length Constraint", GH_Kernel.GH_ParamAccess.item);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialisation ----- //

            Types_Var.Gh_Variable start = null;
            Types_Var.Gh_Variable end = null;
            double length = 0d;

            double weight = 0d;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref start)) { return; };
            if (!DA.GetData(1, ref end)) { return; };
            if (!DA.GetData(2, ref length)) { return; };

            if (!DA.GetData(3, ref weight)) { weight = 1d; };

            // ----- Core ----- //

            /* To Do : Verify that start, end and vector have the same dimentsion. */

            int dimension = start.Value.Dimension;
            GP.Variable[] variables = new GP.Variable[2] { start.Value, end.Value };

            SegmentLength constraintType = new SegmentLength(dimension, length);
            GP.Constraint constraint = new GP.Constraint(constraintType, variables, weight);

            // ----- Set Output ----- //

            DA.SetData(0, constraint);
        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{82CC4A44-328C-4D4A-9C29-8D09290F68D0}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Segment;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }


    /// <summary>
    /// Constraint enforcing a segment to equal a fixed target length. The list of variables for this constraint consists of:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>P<sub>s</sub></term>
    ///         <description> Variable representing the start point of the segment.</description>
    ///     </item>
    ///     <item> 
    ///         <term>P<sub>e</sub></term>
    ///         <description> Variable representing the end point of the segment.</description>
    ///     </item>
    /// </list>
    /// </summary>
    public class SegmentLength : GP.Abstracts.ConstraintType
    {
        /// <summary>
        /// Initialises a new instance of the <see cref="SegmentLength"/> class.
        /// </summary>
        /// <param name="dimension"> Dimension of the variables representing the start and the end segment. </param>
        /// <param name="length"> Value of the target lenght. </param>
        public SegmentLength(int dimension, double length)
        {
            // Start Point : (xs; ys)
            // End Point : (xe, ye)
            // Target Length : Lt
            // (xs - xe)² + (ys - ye)² - Lt² = (xs² - 2 xs xe + xe²) + (ys² - 2 ys ye + ye²) - Lt²

            LinAlg_Mat.Storage.DictionaryOfKeys dok = new LinAlg_Mat.Storage.DictionaryOfKeys(4 * dimension);
            for (int i = 0; i < dimension; i++)
            {
                dok.Add(2d, i, i); dok.Add(-2d, dimension + i, i);
                dok.Add(-2d, i, dimension + i); dok.Add(2d, dimension + i, dimension + i);
            }

            LocalHi = new LinAlg_Mat.Sparse.CompressedColumn(2 * dimension, 2 * dimension, dok);

            LocalBi = null;

            Ci = -(length * length);
        }
    }
}
