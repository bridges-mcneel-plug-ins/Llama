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
    public class Comp_NodeFairness : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="NodeFairness"/> class.
        /// </summary>
        public Comp_NodeFairness()
          : base("Node Fairness", "Fairness",
              "Create a Node Fairness energy for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Energies])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Central Node", "C", "Variable representing the central node to fair. ", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Neighbouring Nodes", "N", "Variable representing the neighouring nodes. ", GH_Kernel.GH_ParamAccess.list);

            pManager.AddNumberParameter("Weight", "W", "Weight of the energy.", GH_Kernel.GH_ParamAccess.item);

            pManager[2].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Energies.Param_Energy(), "Energy", "E", "Node Fairness Energy", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/

            Types.Variables.Gh_Variable node = null;
            List<Types.Variables.Gh_Variable> neighbours =  new List<Types.Variables.Gh_Variable>();


            double weight = 0.0;

            /******************** Get Inputs ********************/

            if (!DA.GetData(0, ref node)) { return; };
            if (!DA.GetDataList(1, neighbours)) { return; };

            if (!DA.GetData(2, ref weight)) { weight = 1.0; };

            /******************** Core ********************/

            int valency = neighbours.Count;
            NodeFairness energyType = new NodeFairness(node.Value.Dimension, valency);

            GP.Variable[] variables = new GP.Variable[valency + 1];
            variables[0] = node.Value;
            for (int i = 0; i < valency; i++) { variables[i + 1] = neighbours[i].Value; }

            GP.Energy energy = new GP.Energy(energyType, variables, weight);
            Types.Energies.Gh_Energy gh_Energy = new Types.Energies.Gh_Energy(energy);

            /******************** Set Output ********************/

            DA.SetData(0, gh_Energy);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{ABD389D2-82B9-4D69-9DEA-B1FF1E86EA5D}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Others;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ComponentAttributes(this);
        }

        #endregion
    }

    // Node car dans le cas ou valence = 4 , on definit deux energies, une dans chaque direction.
    // On est localement dans le lissage de deux polyligne et non d'un maillage


    /// <summary>
    /// Energy enforcing the fairness of a node with regards to its neighbouring nodes. The list of variables of this energy consists in:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>V</term>
    ///         <description> Variable representing the central node.</description>
    ///     </item>
    ///     <item> 
    ///         <term>N<sub>1</sub></term>
    ///         <description> Variables representing the first neighbouring node.</description>
    ///     </item>
    ///     <item> 
    ///         ...
    ///     </item>
    ///     <item> 
    ///         <term>N<sub>n</sub></term>
    ///         <description> Variables representing the last neighbouring node.</description>
    ///     </item>
    /// </list>
    /// </summary>
    /// <remarks> The fairness energy is based on the second order differences of the Laplacian graph. </remarks>
    public class NodeFairness : GP.Abstracts.EnergyType
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="NodeFairness"/> class.
        /// </summary>
        /// <param name="dimension"> Dimension of the node variables. </param>
        /// <param name="valency"> Valency of the node. </param>
        public NodeFairness(int dimension, int valency)
        {
            // ----- Define Ki ----- //

            int count = dimension * (valency + 1);

            int[] rowIndices = new int[count];
            for (int i = 0; i < count; i++)
            {
                rowIndices[i] = i;
            }

            double weight = -(1 / valency);
            double[] values = new double[count];
            for (int i = 0; i < dimension; i++) { values[i] = 1; }
            for (int i = dimension; i < count; i++) { values[i] = weight; }


            LocalKi = new LinAlg_Vect.SparseVector(count, rowIndices, values);

            // ----- Define Si ----- //

            Si = 0d;
        }

        #endregion
    }
}
