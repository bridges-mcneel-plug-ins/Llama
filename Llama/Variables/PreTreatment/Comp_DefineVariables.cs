using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;
using GH_Data = Grasshopper.Kernel.Data;
using GH_Types = Grasshopper.Kernel.Types;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Typ = Llama.Types;
using Params = Llama.Parameters;


namespace Llama.Variables.PreTreatment
{
    /// <summary>
    /// A grasshopper component creating a <see cref="Typ.Variables.Gh_VariableSet"/>.
    /// </summary>
    public class Comp_DefineVariables : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_DefineVariables"/> class.
        /// </summary>
        public Comp_DefineVariables()
          : base("Define Variables", "Set",
              "Defines a set of variables for the Guided Projection algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Variables])
        {
            /* Do Nothing */
        }

        #endregion

        #region Other Methods

        /// <summary>
        /// Creates a list of <see cref="GP.Variable"/> from lists of numerical values.
        /// </summary>
        /// <param name="numbers"> Data tree of numerical values. </param>
        /// <returns> The list of <see cref="GP.Variable"/>. </returns>
        /// <exception cref="ArgumentNullException"> The collection of numerical values cannot be empty. </exception>
        /// <exception cref="ArgumentException"></exception>
        private List<GP.Variable> CreateVariables(GH_Data.GH_Structure<GH_Types.GH_Number> numbers)
        {
            // Verifications
            if (numbers.IsEmpty) { throw new ArgumentNullException(nameof(numbers), "The collection of numerical values cannot be empty."); }


            List<GP.Variable> variables;

            // There are two possible interpretations of the collection of numerical values.
            // - If the input GH_Structure is actually a list (it as only one branch), then each variable corresponds to an item of the list.
            // - If the input GH_Structure is actually a tree (it as more than one branch), the each variable corresponds to a branch of the tree.

            if (numbers.Branches.Count == 1)
            {
                int variableCount = numbers[0].Count;
                variables = new List<GP.Variable>(numbers.DataCount);

                List<GH_Types.GH_Number> list = numbers[0];
                for (int i = 0; i < variableCount; i++)
                {
                    GP.Variable variable = new GP.Variable(list[i].Value);
                    variables.Add(variable);
                }
            }
            else
            {
                int variableCount = numbers.Branches.Count;
                variables = new List<GP.Variable>(variableCount);

                for (int i_Var = 0; i_Var < variableCount; i_Var++)
                {
                    List<GH_Types.GH_Number> branch = numbers[i_Var];

                    double[] components = new double[branch.Count];
                    for (int i_Comp = 0; i_Comp < branch.Count; i_Comp++)
                    {
                        components[i_Comp] = branch[i_Comp].Value;
                    }

                    GP.Variable variable = new GP.Variable(components);
                    variables.Add(variable);
                }
            }

            // ----- Set Output ----- //

            Message = "Generic";
            return variables;
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "N", "Name of the set of variables", GH_Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Numeric Values", "V", "Values for the set variables, represented as lists of numerical values.", GH_Kernel.GH_ParamAccess.tree);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_VariableSet(), "Set of variables", "S", "Set of variables initialised with the specified values.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Variables", "V", "Variables composing the newly constructed Set.", GH_Kernel.GH_ParamAccess.list);
        }


        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Get Inputs ----- //

            string name = "";
            DA.GetData(0, ref name);

            if (!DA.GetDataTree(1, out GH_Data.GH_Structure<GH_Types.GH_Number> numbers)) { return; };
            // GetDataTree does not throw an exception if the conversion failed, but changes a Runtime Message Level.
            // To avoid any NullReferenceException, the level of the Runtime Message is checked.
            if (RuntimeMessageLevel == GH_Kernel.GH_RuntimeMessageLevel.Error)
            {
                throw new InvalidCastException("The input could not be casted to GH_Number.");
            }

            // ----- Core ----- //

            List<GP.Variable> variables = CreateVariables(numbers);

            // ----- Set Output ----- //

            Typ.Variables.Gh_VariableSet gh_Set = new Typ.Variables.Gh_VariableSet(variables, name);

            DA.SetData(0, gh_Set);
            DA.SetDataList(1, variables);

        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{21AC192F-2646-48AF-A13C-E25D4ED73EB0}");

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
