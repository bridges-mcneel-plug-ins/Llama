using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Typ = Llama.Types;
using Params = Llama.Parameters;


namespace Llama.Variables.PostTreatment
{
    /// <summary>
    /// A grasshopper component retrieving a <see cref="Typ.Variables.Gh_VariableSet"/> in a model.
    /// </summary>
    public class Comp_RetrieveVariables : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_RetrieveVariables"/> class.
        /// </summary>
        public Comp_RetrieveVariables()
          : base("Retrieve Variables", "Set",
              "Retrives a set of variables in the model into its numerical values.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Variables])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Set Name", "N", "Name of the set of variables to deconstruct", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params.Models.Param_Model(), "GPA Model", "M", "Assembled Model for the Guided Projection Algorithm.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Variables", "V", "Variables of the set.", GH_Kernel.GH_ParamAccess.list);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialise ----- //

            string name = "";
            Typ.Models.Gh_Model model = new Typ.Models.Gh_Model();

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref name)) { return; } ;
            if (!DA.GetData(1, ref model)) { return; };

            // ----- Core ----- //

            if (model.Sets.TryGetValue(name, out List<GP.Variable> variables))
            {
                DA.SetDataList(0, variables);
            }
            else
            {
                AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Error, "The specified name does not correspond to any variable set i the model.");
                return;
            }

            // ----- Set Output ----- //

            

        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{87C3EAC8-D980-4D05-B54D-3E10FCC628C5}");

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
