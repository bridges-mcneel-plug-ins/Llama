using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using Gh_Types_Euc3D = BRIDGES.McNeel.Grasshopper.Types.Geometry.Euclidean3D;
using Gh_Params_Euc3D = BRIDGES.McNeel.Grasshopper.Parameters.Geometry.Euclidean3D;
using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using GH_Kernel = Grasshopper.Kernel;

using Params = Llama.Parameters;


namespace Llama.Energies.Segment
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
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Energies])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Start Variable", "S", "Variable representing the start of the segment", GH_Kernel.GH_ParamAccess.item);
            pManager.AddParameter(new Params.Variables.Param_Variable(), "End Variable", "E", "Variable representing the end of the segment", GH_Kernel.GH_ParamAccess.item);

            pManager.AddParameter(new Gh_Params_Euc3D.Param_Vector(), "Vector", "V", "Vector to which the segment must be orthogonal.", GH_Kernel.GH_ParamAccess.item);

            pManager.AddNumberParameter("Weight", "W", "Weight of the energy.", GH_Kernel.GH_ParamAccess.item);

            pManager[3].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Energies.Param_Energy(), "Energy", "E", "Segment Orthogonality Energy.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/

            Types.Variables.Gh_Variable start = null;
            Types.Variables.Gh_Variable end = null;

            Gh_Types_Euc3D.Gh_Vector vector = null;

            double weight = 0.0;

            /******************** Get Inputs ********************/

            if (!DA.GetData(0, ref start)) { return; };
            if (!DA.GetData(1, ref end)) { return; };

            if (!DA.GetData(2, ref vector)) { return; };

            if (!DA.GetData(3, ref weight)) { weight = 1.0; };

            /******************** Core ********************/

            int dimension = vector.Value.Dimension;
            if (dimension != start.Value.Dimension || dimension != end.Value.Dimension)
            {
                throw new ArgumentException("The start and end variables must have the same number of components than the vector.", new RankException());
            }

            double[] components = new double[] { vector.Value.X, vector.Value.Y, vector.Value.Z };
            GP.EnergyTypes.SegmentOrthogonality energyType = new GP.EnergyTypes.SegmentOrthogonality(components);

            GP.Variable[] variables = new GP.Variable[2] { start.Value, end.Value };

            GP.Energy energy = new GP.Energy(energyType, variables, weight);
            Types.Energies.Gh_Energy gh_Energy = new Types.Energies.Gh_Energy(energy);

            /******************** Set Output ********************/

            DA.SetData(0, gh_Energy);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{9FC632E5-85B5-4857-84D8-B1B22FE7D7FC}");

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
}
