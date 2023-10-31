using System;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Types_Var = Llama.Types.Variables;
using Types_Cstr = Llama.Types.Constraints;
using Params_Var = Llama.Parameters.Variables;
using Params_Cstr = Llama.Parameters.Constraints;


namespace Llama.Constraints.Numeric
{
    /// <summary>
    /// A grasshopper component creating a <see cref="GP.QuadraticConstraintTypes.LowerBound"/>-based <see cref="GP.Constraint"/>..
    /// </summary>
    public class Comp_LowerBound : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_LowerBound"/> class.
        /// </summary>
        public Comp_LowerBound()
          : base("Lower Bound", "Lower Bound",
              "Create an Lower Bound constraint for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Constraints])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params_Var.Param_Variable(), "Scalar Variable", "N", "Scalar variable whose value should be larger than a specified value.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Lower Bound", "B", "Lower Bound of the scalar variable.", GH_Kernel.GH_ParamAccess.item);

            pManager.AddNumberParameter("Weight", "W", "Weight of the constraint.", GH_Kernel.GH_ParamAccess.item);

            pManager[2].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params_Cstr.Param_Constraint(), "Constraint", "C", "Vector Length Constraint", GH_Kernel.GH_ParamAccess.item);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialisation ----- //

            Types_Var.Gh_Variable scalar = null;
            double bound = 0d;

            double weight = 0d;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref scalar)) { return; };
            if (!DA.GetData(1, ref bound)) { return; };

            if (!DA.GetData(2, ref weight)) { weight = 1d; };

            // ----- Core ----- //

            /* To Do : Verify that start, end and vector have the same dimentsion. */

            GP.Variable dummy = new GP.Variable(0d);
            GP.Variable[] variables = new GP.Variable[2] { scalar.Value, dummy };

            GP.QuadraticConstraintTypes.LowerBound constraintType = new GP.QuadraticConstraintTypes.LowerBound(bound);
            GP.Constraint constraint = new GP.Constraint(constraintType, variables, weight);

            // ----- Set Output ----- //

            DA.SetData(0, constraint);
        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{B885A300-C37E-41AF-A441-47F5FC25E8B7}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Numeric;

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
