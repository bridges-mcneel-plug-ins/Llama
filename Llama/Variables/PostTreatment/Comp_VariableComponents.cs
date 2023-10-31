using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Params = Llama.Parameters;


namespace Llama.Variables.PostTreatment
{
    /// <summary>
    /// A grasshopper component translating <see cref="GP.Variable"/> to their numerical components.
    /// </summary>
    public class Comp_VariableComponents : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_VariableComponents"/> class.
        /// </summary>
        public Comp_VariableComponents()
          : base("Variable Components", "Comp.",
              "Translates a variable into its numerical components.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Variables])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Variable", "V", "Variable to translate into its components.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Components", "C", "Components of the variable.", GH_Kernel.GH_ParamAccess.list);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            GP.Variable variable = null;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref variable)) { return; }

            // ----- Core ----- //
            List<double> components = new List<double>(variable.Dimension);
            for (int i = 0; i < variable.Dimension; i++)
            {
                components.Add(variable[i]);
            }

            // ----- Set Output ----- //

            DA.SetDataList(0, components);

        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{5440C504-9B4F-4CFF-83D9-6B325E74BC8A}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return null;
            }
        }

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.PostTreatment;


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }
}
