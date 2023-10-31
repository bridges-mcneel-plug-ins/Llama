using System;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Typ = Llama.Types.Models;


namespace Llama.Parameters.Models
{
    /// <summary>
    /// A <see cref="Typ.Gh_Model"/> grasshopper parameter.
    /// </summary>
    public class Param_Model : GH_Kernel.GH_Param<Typ.Gh_Model>
    {
        #region Properties

        /// <inheritdoc cref="GH_Kernel.IGH_PreviewObject.Hidden"/>
        public bool Hidden => true;

        /// <inheritdoc cref="GH_Kernel.IGH_PreviewObject.IsPreviewCapable"/>
        public bool IsPreviewCapable => false;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new instance of <see cref="Param_Model"/>.
        /// </summary>
        public Param_Model()
          : base("Guided Projection Algorithm Model", "GPA Model", "Contains a collection of Guided Projection Algorithm model.",
                Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Parameters], GH_Kernel.GH_ParamAccess.item)
        {
            /* Do Nothing */

        }

        #endregion


        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{309F1C6F-0D58-4D39-A680-810037853971}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Models;

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
