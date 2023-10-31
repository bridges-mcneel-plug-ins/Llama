using System;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Types_Var = Llama.Types.Variables;
using Types_Cstr = Llama.Types.Constraints;
using Params_Var = Llama.Parameters.Variables;
using Params_Cstr = Llama.Parameters.Constraints;


namespace Llama.Constraints.Vector
{
    /// <summary>
    /// A grasshopper component creating a <see cref="GP.QuadraticConstraintTypes.VectorLength"/>-based <see cref="GP.Constraint"/>..
    /// </summary>
    public class Comp_VectorLength : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_VectorLength"/> class.
        /// </summary>
        public Comp_VectorLength()
          : base("Vector Length", "Length",
              "Create a Vector Length constraint for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Constraints])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params_Var.Param_Variable(), "Vector Variable", "V", "Vector-like variable whose length should equal a specified target length", GH_Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Taget Length", "L", "Target length of the vector variable.", GH_Kernel.GH_ParamAccess.item);

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

            Types_Var.Gh_Variable vector = null;
            double length = 0d;

            double weight = 0d;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref vector)) { return; };
            if (!DA.GetData(1, ref length)) { return; };

            if (!DA.GetData(2, ref weight)) { weight = 1d; };

            // ----- Core ----- //

            /* To Do : Verify that start, end and vector have the same dimentsion. */

            int dimension = vector.Value.Dimension;
            GP.Variable[] variables = new GP.Variable[1] { vector.Value };

            GP.QuadraticConstraintTypes.VectorLength constraintType = new GP.QuadraticConstraintTypes.VectorLength(length, dimension);
            GP.Constraint constraint = new GP.Constraint(constraintType, variables, weight);

            // ----- Set Output ----- //

            DA.SetData(0, constraint);
        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{3925C384-0698-405D-9C77-ACD0A85AC720}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Vector;

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
