using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Types = Llama.Types;
using Params = Llama.Parameters;


namespace Llama.Models
{
    /// <summary>
    /// A grasshopper component creating a <see cref="GP.GuidedProjectionAlgorithm"/> model.
    /// </summary>
    public class Comp_AssembleModel : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_AssembleModel"/> class.
        /// </summary>
        public Comp_AssembleModel()
          : base("Assemble Guided Projection Model", "Assemble GP",
              "Create a model for the Guided Projection Algorithm (GPA).",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Solvers])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_VariableSet(), "Sets", "S", "Sets of variable used by the energies and constraints.", GH_Kernel.GH_ParamAccess.list);
            pManager.AddParameter(new Params.Energies.Param_Energy(), "Energies", "E", "Energies", GH_Kernel.GH_ParamAccess.list);
            pManager.AddParameter(new Params.Constraints.Param_Constraint(), "Constraints", "C", "Constraints", GH_Kernel.GH_ParamAccess.list);

            pManager[1].Optional = true;
            pManager[2].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Models.Param_Model(), "Model", "M", "Model for the Guided Projection Algorithm.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ---------- Initialisation ---------- //

            List<Types.Variables.Gh_VariableSet> gh_Sets = new List<Types.Variables.Gh_VariableSet>();
            List<Types.Energies.Gh_Energy> gh_Energies = new List<Types.Energies.Gh_Energy>();
            List<Types.Constraints.Gh_Constraint> gh_Constraints = new List<Types.Constraints.Gh_Constraint>();

            // ---------- Get Inputs ---------- //

            if (!DA.GetDataList(0, gh_Sets)) { return; };

            bool hasEnergies = DA.GetDataList(1, gh_Energies);
            bool hasConstraints = DA.GetDataList(2, gh_Constraints);
            if ((!hasEnergies || gh_Energies.Count == 0) & (!hasConstraints || gh_Constraints.Count == 0))
            {
                AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Warning, "Input parameters E and C failed to collect data. The model needs energies or constraint to run.");


                Params.Output[0].Phase = GH_Kernel.GH_SolutionPhase.Blank;
                Params.Output[0].ClearData();


                // this.Params.Output[0].ExpireSolution(); ????

                return;
            }

            // ---------- Core ---------- //

            int setVariableCount = 0;
            foreach (Types.Variables.Gh_VariableSet gh_Set in gh_Sets)
            {
                setVariableCount += gh_Set.Count;
            }

            GP.GuidedProjectionAlgorithm gpa = new GP.GuidedProjectionAlgorithm(1e-4, 100);


            // ---------- Variables ---------- //


            //
            //  Il faut copier les variables pour de problème de reference a priori ...
            //

            // Key : Old Variables
            // Value : Copied Variables
            Dictionary<GP.Variable, GP.Variable> oldToNew = new Dictionary<GP.Variable, GP.Variable>(setVariableCount);


            Dictionary<string, List<GP.Variable>> sets = new Dictionary<string, List<GP.Variable>>(gh_Sets.Count);


            foreach (Types.Variables.Gh_VariableSet gh_Set in gh_Sets)
            {
                List<GP.Variable> variables = new List<GP.Variable>(gh_Set.Count);
                foreach (GP.Variable variable in gh_Set)
                {
                    if (oldToNew.ContainsKey(variable))
                    {
                        AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Warning, $"The set {gh_Set.Name} contains a variable that was already added to the model.");
                    }
                    else
                    {
                        double[] components = variable.ToArray();
                        GP.Variable newVariable = gpa.AddVariable(components);

                        oldToNew.Add(variable, newVariable);

                        variables.Add(newVariable);
                    }
                }

                sets.Add(gh_Set.Name, variables);
            }

            // ---------- Energies ---------- //

            foreach (Types.Energies.Gh_Energy gh_Energy in gh_Energies)
            {
                IReadOnlyList<GP.Variable> variables = gh_Energy.Value.Variables;

                List<GP.Variable> newVariables = new List<GP.Variable>(variables.Count);

                foreach (GP.Variable variable in variables)
                {
                    if (oldToNew.ContainsKey(variable))
                    {
                        GP.Variable newVariable = oldToNew[variable];
                        newVariables.Add(newVariable);
                    }
                    else
                    {
                        AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Remark, $"A variable which does not belong to any set was added to the model. It might be a dummy variable.");

                        double[] components = variable.ToArray();
                        GP.Variable newVariable = gpa.AddVariable(components);

                        oldToNew.Add(variable, newVariable);

                        newVariables.Add(newVariable);
                    }
                }

                gpa.AddEnergy(gh_Energy.Value.Type, newVariables, gh_Energy.Value.Weight);

                // if (!isEnergyAdded)
                // {
                //     this.AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Warning, $"This energy is a duplicate, it was already added to the model.");
                // }
            }

            // ---------- Constraints ---------- //

            foreach (Types.Constraints.Gh_Constraint gh_Constraint in gh_Constraints)
            {
                IReadOnlyList<GP.Variable> variables = gh_Constraint.Value.Variables;

                List<GP.Variable> newVariables = new List<GP.Variable>(variables.Count);

                foreach (GP.Variable variable in variables)
                {
                    if (oldToNew.ContainsKey(variable))
                    {
                        GP.Variable newVariable = oldToNew[variable];
                        newVariables.Add(newVariable);
                    }
                    else
                    {
                        AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Remark, $"A variable which does not belong to any set was added to the model. It might be a dummy variable.");

                        double[] components = variable.ToArray();
                        GP.Variable newVariable = gpa.AddVariable(components);

                        oldToNew.Add(variable, newVariable);

                        newVariables.Add(newVariable);
                    }
                }

                gpa.AddConstraint(gh_Constraint.Value.Type, newVariables, gh_Constraint.Value.Weight);

                // if (!isConstraintAdded)
                // {
                //     this.AddRuntimeMessage(GH_Kernel.GH_RuntimeMessageLevel.Warning, $"This constraint is a duplicate, it was already added to the model.");
                // }
            }


            // ---------- Set Output ---------- //

            Types.Models.Gh_Model gh_Model = new Types.Models.Gh_Model(gpa, ref sets);

            DA.SetData(0, gh_Model);
        }

        #endregion

        #region Override : GH_DocumentObject

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{0DB58DA4-D469-4AF7-AC13-D9B630127B42}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Solver;


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
