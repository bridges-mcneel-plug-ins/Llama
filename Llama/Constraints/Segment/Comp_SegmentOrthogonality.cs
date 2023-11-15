using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;
using LinAlg_Mat = BRIDGES.LinearAlgebra.Matrices;

using Gh_Types_Euc3D = BRIDGES.McNeel.Grasshopper.Types.Geometry.Euclidean3D;
using Gh_Params_Euc3D = BRIDGES.McNeel.Grasshopper.Parameters.Geometry.Euclidean3D;
using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using GH_Kernel = Grasshopper.Kernel;

using Params = Llama.Parameters;
using Rhino.Geometry;


namespace Llama.Constraints.Segment
{
    /// <summary>
    /// A grasshopper component creating a <see cref="GP.EnergyTypes.SegmentOrthogonality"/>-based <see cref="GP.Energy"/>.
    /// </summary>
    public class Comp_SegmentOrthogonality : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_SegmentOrthogonality"/> class.
        /// </summary>
        public Comp_SegmentOrthogonality()
          : base("Segment Orthogonality", "Ortho.",
              "Create a Segment Orthogonality energy for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Constraints])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Start Variable", "S", "Variable representing the start of the segment.", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params.Variables.Param_Variable(), "End Variable", "E", "Variable representing the end of the segment.", GH_Kernel.GH_ParamAccess.item);

            pManager.AddParameter(new Params.Variables.Param_Variable(), "Vector", "V", "Variable representing the vector to which the segment must be orthogonal.", GH_Kernel.GH_ParamAccess.item);


            pManager.AddNumberParameter("Weight", "W", "Weight of the constraint.", GH_Kernel.GH_ParamAccess.item);

            pManager[3].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Constraints.Param_Constraint(), "Constraint", "C", "Segment Orthogonality Constraint.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/

            Types.Variables.Gh_Variable start = null;
            Types.Variables.Gh_Variable end = null;
            Types.Variables.Gh_Variable normal = null;


            double weight = 0.0;

            /******************** Get Inputs ********************/

            if (!DA.GetData(0, ref start)) { return; };
            if (!DA.GetData(1, ref end)) { return; };

            if (!DA.GetData(2, ref normal)) { return; };

            if (!DA.GetData(3, ref weight)) { weight = 1.0; };

            /******************** Core ********************/

            int dimension = normal.Value.Dimension;
            if (dimension != start.Value.Dimension || dimension != end.Value.Dimension)
            {
                throw new ArgumentException("The start and end variables must have the same number of components than the vector variable.", new RankException());
            }

            SegmentOrthogonality constraintType = new SegmentOrthogonality(dimension);

            GP.Variable[] variables = new GP.Variable[3] { start.Value, end.Value, normal.Value };

            GP.Constraint constraint = new GP.Constraint(constraintType, variables, weight);
            Types.Constraints.Gh_Constraint gh_Constraint = new Types.Constraints.Gh_Constraint(constraint);

            /******************** Set Output ********************/

            DA.SetData(0, gh_Constraint);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{EEFF6016-F2B7-4A66-B798-D419CB3F97F4}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Segment;

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
    /// Constraint enforcing a segment to be orthogonal to a variable direction. The list of variables for this constraint consists of:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>P<sub>s</sub></term>
    ///         <description> Variable representing the start point of the segment.</description>
    ///     </item>
    ///     <item> 
    ///         <term>P<sub>e</sub></term>
    ///         <description> Variable representing the end point of the segment.</description>
    ///     </item>
    ///     <item> 
    ///         <term>V</term>
    ///         <description> Variable representing the direction to which the segment must be orthogonal.</description>
    ///     </item>
    /// </list>
    /// </summary>
    public class SegmentOrthogonality : GP.Abstracts.ConstraintType
    {
        /// <summary>
        /// Initialises a new instance of the <see cref="SegmentLength"/> class.
        /// </summary>
        /// <param name="dimension"> Dimension of the variables representing the start and the end segment. </param>
        public SegmentOrthogonality(int dimension)
        {

            LinAlg_Mat.Storage.DictionaryOfKeys dok = new LinAlg_Mat.Storage.DictionaryOfKeys(4 * dimension);
            for (int i = 0; i < dimension; i++)
            {
                dok.Add(-1d, (2 * dimension) + i, i); dok.Add(-1d, i, (2 * dimension) + i);
                dok.Add(1d, (2 * dimension) + i, dimension + i); dok.Add(1d, dimension + i, (2 * dimension) + i);
            }

            LocalHi = new LinAlg_Mat.Sparse.CompressedColumn(3 * dimension, 3 * dimension, dok);

            LocalBi = null;

            Ci = 0;
        }
    }

}
