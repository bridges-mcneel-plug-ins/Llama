using System;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Types_Var = Llama.Types.Variables;
using Types_Cstr = Llama.Types.Constraints;
using Params_Var = Llama.Parameters.Variables;
using Params_Cstr = Llama.Parameters.Constraints;


namespace Llama.Constraints
{
    /// <summary>
    /// A grasshopper component creating a <see cref="GP.QuadraticConstraintTypes.CoherentLength"/>-based <see cref="GP.Constraint"/>..
    /// </summary>
    public class Comp_CoherentLength : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_CoherentLength"/> class.
        /// </summary>
        public Comp_CoherentLength()
          : base("Coherent Length", "Coh. Length",
              "Create a Coherent Length constraint for the Guided Projection Algorithm.",
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
            pManager.AddParameter(new Params_Var.Param_Variable(), "Length Variable", "L", "Variable representing the distance between the start and the end of the segment", GH_Kernel.GH_ParamAccess.item);

            pManager.AddNumberParameter("Weight", "W", "Weight of the constraint.", GH_Kernel.GH_ParamAccess.item);

            pManager[3].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params_Cstr.Param_Constraint(), "Constraint", "C", "Coherent Length Constraint", GH_Kernel.GH_ParamAccess.item);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialisation ----- //

            Types_Var.Gh_Variable start = null;
            Types_Var.Gh_Variable end = null;
            Types_Var.Gh_Variable length = null;

            double weight = 0d;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref start)) { return; };
            if (!DA.GetData(1, ref end)) { return; };
            if (!DA.GetData(2, ref length)) { return; };

            if (!DA.GetData(3, ref weight)) { weight = 1d; };

            // ----- Core ----- //

            /* To Do : Verify that start, end and vector have the same dimentsion. */

            int dimension = start.Value.Dimension;
            GP.Variable[] variables = new GP.Variable[3] { start.Value, end.Value, length.Value };

            GP.QuadraticConstraintTypes.CoherentLength constraintType = new GP.QuadraticConstraintTypes.CoherentLength(dimension);
            GP.Constraint constraint = new GP.Constraint(constraintType, variables, weight);

            // ----- Set Output ----- //

            DA.SetData(0, constraint);
        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{DD82A6D7-ABFA-48A1-9A86-858015DA2A2D}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Others;

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
}
