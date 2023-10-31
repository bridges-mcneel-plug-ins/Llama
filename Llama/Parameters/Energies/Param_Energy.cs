using System;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Typ = Llama.Types.Energies;


namespace Llama.Parameters.Energies
{
    /// <summary>
    /// A <see cref="Typ.Gh_Energy"/> grasshopper parameter.
    /// </summary>
    public class Param_Energy : GH_Kernel.GH_Param<Typ.Gh_Energy>
    {
        #region Properties

        /// <inheritdoc cref="GH_Kernel.IGH_PreviewObject.Hidden"/>
        public bool Hidden => true;

        /// <inheritdoc cref="GH_Kernel.IGH_PreviewObject.IsPreviewCapable"/>
        public bool IsPreviewCapable => false;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new instance of <see cref="Param_Energy"/>.
        /// </summary>
        public Param_Energy()
          : base("Energy", "Energy", "Contains a collection of energies for the Guided Projection Algorithm.",
                Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Parameters], GH_Kernel.GH_ParamAccess.item)
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{D084ADB9-E601-4208-97DB-DE57063A7A0C}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Energies;

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Icon"/>
        protected override System.Drawing.Bitmap Icon => null;


        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.CreateAttributes()"/>
        public override void CreateAttributes()
        {
            m_attributes = new Gh_Disp_Euc3D.ParameterAttributes(this);
        }

        #endregion

        #region Override : GH_ActiveObject

        /// <inheritdoc cref="GH_Kernel.GH_ActiveObject.Locked"/>
        public override bool Locked
        {
            get { return base.Locked; }
            set
            {
                if (base.Locked != value)
                {
                    base.Locked = value;
                }
            }
        }

        #endregion
    }
}
