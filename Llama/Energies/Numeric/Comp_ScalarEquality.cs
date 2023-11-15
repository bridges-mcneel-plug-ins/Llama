using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;
using LinAlg_Vect = BRIDGES.LinearAlgebra.Vectors;

using Gh_Types_Euc3D = BRIDGES.McNeel.Grasshopper.Types.Geometry.Euclidean3D;
using Gh_Params_Euc3D = BRIDGES.McNeel.Grasshopper.Parameters.Geometry.Euclidean3D;
using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using GH_Kernel = Grasshopper.Kernel;

using Params = Llama.Parameters;


namespace Llama.Energies
{
    /// <summary>
    /// A grasshopper component creating a <see cref="NodeFairness"/>-based <see cref="GP.Energy"/>.
    /// </summary>
    public class Comp_ScalarEquality : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="NodeFairness"/> class.
        /// </summary>
        public Comp_ScalarEquality()
          : base("Scalar Equality", "Equality",
              "Create a Scalar Equality energy for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Energies])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Scalar Variable", "V", "Variable representing the scalar value. ", GH_Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("Target Value", "T", "Target Value for the scalar variable.", GH_Kernel.GH_ParamAccess.item);

            pManager.AddNumberParameter("Weight", "W", "Weight of the energy.", GH_Kernel.GH_ParamAccess.item);

            pManager[2].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Energies.Param_Energy(), "Energy", "E", "Scalar Equality Energy", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            // ----- Initialisation ----- //

            Types.Variables.Gh_Variable variable = null;
            double targetValue = 0d;

            double weight = 0d;

            // ----- Get Inputs ----- //

            if (!DA.GetData(0, ref variable)) { return; };
            if (!DA.GetData(1, ref targetValue)) { return; };

            if (!DA.GetData(2, ref weight)) { weight = 1.0; };

            // ----- Core ----- //

            ScalarEquality energyType = new ScalarEquality(targetValue);

            GP.Variable[] variables = new GP.Variable[1] { variable.Value };

            GP.Energy energy = new GP.Energy(energyType, variables, weight);
            Types.Energies.Gh_Energy gh_Energy = new Types.Energies.Gh_Energy(energy);

            // ----- Set Output ----- //

            DA.SetData(0, gh_Energy);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{7D10E61F-BD86-4DB7-A286-9266EE723ED5}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Numeric;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

    /// <summary>
    /// Energy enforcing a numeric variable to be equal to a fixed value. The list of variables of this energy consists in:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>V</term>
    ///         <description> Scalar variable. which should match the specified fixed value. </description>
    ///     </item>
    /// </list>
    /// </summary>
    public class ScalarEquality : GP.Abstracts.EnergyType
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="ScalarEquality"/> class.
        /// </summary>
        /// <param name="value"> Target value for the numeric variable. </param>
        public ScalarEquality(double value)
        {
            // ----- Define Ki ----- //

            int[] rowIndices = new int[1] { 0 };
            double[] values = new double[1] { 1d };

            LocalKi = new LinAlg_Vect.SparseVector(1, rowIndices, values);

            // ----- Define Si ----- //

            Si = value;
        }

        #endregion
    }
}
