using System;
using System.Collections.Generic;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Types = Grasshopper.Kernel.Types;


namespace Llama.Types.Models
{
    /// <summary>
    /// Class defining a grasshopper type for an <see cref="GP.GuidedProjectionAlgorithm"/>.
    /// </summary>
    public class Gh_Model : GH_Types.GH_Goo<GP.GuidedProjectionAlgorithm>
    {
        #region Fields

        /// <summary>
        /// Dictionnary containing the set of variables, identified by their name.
        /// </summary>
        public Dictionary<string, List<GP.Variable>> Sets;

        #endregion


        #region Constructors

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Model" /> class.
        /// </summary>
        public Gh_Model() { /* Do Nothing */ }

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Model" /> class from another <see cref="Gh_Model"/>..
        /// </summary>
        /// <param name="gh_Model"> <see cref="Gh_Model"/> to duplicate. </param>
        public Gh_Model(Gh_Model gh_Model)
        {
            Value = gh_Model.Value;
        }

        /// <summary>
        /// Initialises a new instance of <see cref="Gh_Model"/> class from a <see cref="GP.GuidedProjectionAlgorithm"/>.
        /// </summary>
        /// <param name="gpa"> <see cref="GP.GuidedProjectionAlgorithm"/> for the grasshopper type.</param>
        public Gh_Model(GP.GuidedProjectionAlgorithm gpa)
        {
            Value = gpa;
        }

        /// <summary>
        /// Initialises a new instance of <see cref="Gh_Model"/> class from a <see cref="GP.GuidedProjectionAlgorithm"/>.
        /// </summary>
        /// <param name="gpa"> <see cref="GP.GuidedProjectionAlgorithm"/> for the grasshopper type.</param>
        /// <param name="sets"> Sets of variables used in the Guided Projection Algorithm. </param>
        public Gh_Model(GP.GuidedProjectionAlgorithm gpa, ref Dictionary<string, List<GP.Variable>> sets)
        {
            Value = gpa;
            Sets = sets;
        }

        #endregion


        #region Override : GH_Goo<>

        // ---------- Properties ----------//

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.IsValid"/>
        public override bool IsValid => !(Value is null);

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.TypeDescription"/>
        public override string TypeDescription => $"Grasshopper type containing a {typeof(GP.GuidedProjectionAlgorithm)}.";

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.TypeName"/>
        public override string TypeName => nameof(Gh_Model);


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.Duplicate"/>
        public override GH_Types.IGH_Goo Duplicate() => new Gh_Model(this);


        /// <inheritdoc cref="GH_Types.GH_Goo{T}.CastFrom(object)"/>
        public override bool CastFrom(object source)
        {
            if (source == null) { return false; }

            var type = source.GetType();


            // ----- Otherwise ----- //

            return false;
        }

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.CastTo{Q}(ref Q)"/>
        public override bool CastTo<T>(ref T target)
        {
            // ----- Otherwise ----- //

            return false;
        }

        #endregion

        #region Override : Object

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.ToString"/>
        public override string ToString()
        {
            return Value is null ? "Empty Gh_Model parameter" : $"(V:{Value.ComponentCount}, E:{Value.EnergyCount}, C:{Value.ConstraintCount})";
        }

        #endregion
    }
}
