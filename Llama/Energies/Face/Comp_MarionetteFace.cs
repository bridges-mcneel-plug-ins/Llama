using System;
using System.Collections.Generic;

using Euc3D = BRIDGES.Geometry.Euclidean3D;
using GP = BRIDGES.Solvers.GuidedProjection;
using LinAlg_Vect = BRIDGES.LinearAlgebra.Vectors;

using Gh_Types_Euc3D = BRIDGES.McNeel.Grasshopper.Types.Geometry.Euclidean3D;
using Gh_Params_Euc3D = BRIDGES.McNeel.Grasshopper.Parameters.Geometry.Euclidean3D;
using Gh_Disp_Euc3D = BRIDGES.McNeel.Grasshopper.Display.Geometry.Euclidean3D;

using GH_Kernel = Grasshopper.Kernel;

using Params = Llama.Parameters;


namespace Llama.Energies.Face
{
    /// <summary>
    /// A grasshopper component creating a <see cref="GP.EnergyTypes.SegmentParallelity"/>-based <see cref="GP.Energy"/>.
    /// </summary>
    public class Comp_MarionetteFace : GH_Kernel.GH_Component
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="Comp_MarionetteFace"/> class.
        /// </summary>
        public Comp_MarionetteFace()
          : base("Marionette Face", "Marionette",
              "Create a Marionette Face energy for the Guided Projection Algorithm.",
              Settings.CategoryName, Settings.SubCategoryName[Llama.SubCategory.Energies])
        {
            /* Do Nothing */
        }

        #endregion


        #region Override : GH_Component

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterInputParams(GH_InputParamManager)"/>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddParameter(new Gh_Params_Euc3D.Param_Point(), "Points", "P", "Points of the face", GH_Kernel.GH_ParamAccess.list);
            pManager.AddParameter(new Params.Variables.Param_Variable(), "Height Variables", "H", "Variables representing the heights of the face points.", GH_Kernel.GH_ParamAccess.list);

            pManager.AddNumberParameter("Weight", "W", "Weight of the energy.", GH_Kernel.GH_ParamAccess.item);

            pManager[2].Optional = true;
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.RegisterOutputParams(GH_OutputParamManager)"/>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddParameter(new Params.Energies.Param_Energy(), "Energy", "E", "Face Marionette Energy.", GH_Kernel.GH_ParamAccess.item);
        }

        /// <inheritdoc cref="GH_Kernel.GH_Component.SolveInstance(GH_Kernel.IGH_DataAccess)"/>
        protected override void SolveInstance(GH_Kernel.IGH_DataAccess DA)
        {
            /******************** Initialisation ********************/
            List<Euc3D.Point> points = new List<Euc3D.Point>(4);
            List<GP.Variable> variables = new List<GP.Variable>(4);

            double weight = 0d;

            /******************** Get Inputs ********************/

            if (!DA.GetData(0, ref points)) { return; };
            if (!DA.GetData(1, ref variables)) { return; };

            if (!DA.GetData(2, ref weight)) { weight = 1d; };

            /******************** Core ********************/

            if (points.Count != 4 || variables.Count != 4)
            {
                throw new ArgumentException("Exactly four points and four variables must be provided.", new RankException());
            }

            MarionnetteFace energyType = new MarionnetteFace(points[0], points[1], points[2], points[3]);

            GP.Energy energy = new GP.Energy(energyType, variables, weight);
            Types.Energies.Gh_Energy gh_Energy = new Types.Energies.Gh_Energy(energy);

            /******************** Set Output ********************/

            DA.SetData(0, gh_Energy);
        }

        #endregion

        #region Override : GH_DocumentObject

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.ComponentGuid"/>
        public override Guid ComponentGuid => new Guid("{A29181AD-94F6-4C8B-9C64-075795C8B994}");

        /// <inheritdoc cref="GH_Kernel.GH_DocumentObject.Exposure"/>
        public override GH_Kernel.GH_Exposure Exposure => (GH_Kernel.GH_Exposure)TabExposure.Face;

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
    /// Energy enforcing quadrilateral face planarity on the height field of the vertices. The list of variables of this energy consists in:
    /// <list type="bullet">
    ///     <item> 
    ///         <term>z<sub>a</sub></term>
    ///         <description> Height coordinate of the first point. </description>
    ///     </item>
    ///     <item> 
    ///         <term>z<sub>b</sub></term>
    ///         <description> Height coordinate of the second point. </description>
    ///     </item>
    ///     <item> 
    ///         <term>z<sub>c</sub></term>
    ///         <description> Height coordinate of the third point. </description>
    ///     </item>
    ///     <item> 
    ///         <term>z<sub>d</sub></term>
    ///         <description> Height coordinate of the quadarilateral point. </description>
    ///     </item>
    /// </list>
    /// </summary>
    public class MarionnetteFace : GP.Abstracts.EnergyType
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of the <see cref="MarionnetteFace"/> class.
        /// </summary>
        /// <param name="pointA"> First point of the quadrilateral face. </param>
        /// <param name="pointB"> Second point of the quadrilateral face. </param>
        /// <param name="pointC"> Third point of the quadrilateral face. </param>
        /// <param name="pointD"> Fourth point of the quadrilateral face. </param>
        public MarionnetteFace(Euc3D.Point pointA, Euc3D.Point pointB, Euc3D.Point pointC, Euc3D.Point pointD)
        {
            // ----- Define Ki ----- //

            double d_BC = ((pointB.X - pointA.X) * (pointC.Y - pointA.Y)) - ((pointB.Y - pointA.Y) * (pointC.X - pointA.X));
            double d_BD = ((pointB.X - pointA.X) * (pointD.Y - pointA.Y)) - ((pointB.Y - pointA.Y) * (pointD.X - pointA.X));
            double d_CD = ((pointC.X - pointA.X) * (pointD.Y - pointA.Y)) - ((pointC.Y - pointA.Y) * (pointD.X - pointA.X));

            int[] rowIndices = new int[4] { 0, 1, 2, 3};
            double[] values = new double[4] { -(d_BC + d_BD + d_CD), d_CD, d_BD, d_BC };

            LocalKi = new LinAlg_Vect.SparseVector(4, rowIndices, values);

            // ----- Define Si ----- //

            Si = 0d;
        }

        #endregion
    }
}
