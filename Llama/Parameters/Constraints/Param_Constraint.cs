using System;

using GH_Kernel = Grasshopper.Kernel;

using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using Typ = Llama.Types.Constraints;


namespace Llama.Parameters.Constraints
{
    /// <summary>
    /// A <see cref="Typ.Gh_Constraint"/> grasshopper parameter.
    /// </summary>
    public class Param_Constraint : GH_Kernel.GH_Param<Typ.Gh_Constraint>
    {
        #region Properties

        /// <inheritdoc cref="GH_Kernel.IGH_PreviewObject.Hidden"/>
        public bool Hidden => true;

        /// <inheritdoc cref="GH_Kernel.IGH_PreviewObject.IsPreviewCapable"/>
        public bool IsPreviewCapable => false;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new instance of <see cref="Param_Constraint"/>.
        /// </summary>
        public Param_Constraint()
          : base("Quadratic Constraint", "Q. Constraint", "Contains a collection of quadratic constraints for the Guided Projection Algorithm.",
                Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Parameters], GH_Kernel.GH_ParamAccess.item)
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{2259E43C-6756-4258-8130-DEAB8AACF5EB}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Constraints;

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