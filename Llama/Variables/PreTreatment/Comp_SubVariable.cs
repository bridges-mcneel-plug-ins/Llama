using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using RH_Geo = Rhino.Geometry;

using GH_Kernel = Grasshopper.Kernel;
using GH_Data = Grasshopper.Kernel.Data;
using GH_Types = Grasshopper.Kernel.Types;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Typ = Llama.Types;
using Params = Llama.Parameters;


namespace Llama.Variables
{
    /// <summary>
    /// A grasshopper component creating a subvariable from a <see cref="GP.Variable"/>.
    /// </summary>
    public class Comp_SubVariable : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_SubVariable"/> class.
        /// </summary>
        public Comp_SubVariable()
          : base("Sub-Variables", "Sub V",
              "Defines a SubVariable, i.e. a subset of variable's components.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Variables])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Variables", "V", "Variables whose components are used to create a Sub-Variable.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Start Index", "S", "Zero-based index of the first component's.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddIntegerParameter("Count", "C", "Number of components of the Sub-Variable.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Sub-Variable", "S", "Sub-Variable created from the given components.", GH_Kernel.GH_ParamAccess.item);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Get Inputs ----- //

            GP.Variable variable = null;
            int start = 0, count = 0;

            if (!DA.GetData(0, ref variable)) { return; }
            if (!DA.GetData(1, ref start)) { return; }
            if (!DA.GetData(2, ref count)) { return; }

            // ----- Core ----- //

            if (variable.Dimension < start + count) { throw new ArgumentOutOfRangeException("The given subset is outside the range of the variables components."); }

            // GP.Variable variable = ;

            // ----- Set Output ----- //

            DA.SetData(0, variable);

        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{A54EC1AF-4379-4ED2-997F-89E3018B342A}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return null;
            }
        }

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.PreTreatment;


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }
}
