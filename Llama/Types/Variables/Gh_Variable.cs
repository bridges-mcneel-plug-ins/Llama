using System;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Types = Grasshopper.Kernel.Types;


namespace Llama.Types.Variables
{
    /// <summary>
    /// Class defining a grasshopper type, which represents a variable.
    /// </summary>
    public class Gh_Variable : GH_Types.GH_Goo<GP.Variable>
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Variable" /> class.
        /// </summary>
        public Gh_Variable() { /* Do Nothing */ }

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Variable" /> class from another <see cref="Gh_Variable"/>.
        /// </summary>
        /// <param name="gh_Variable"> <see cref="Gh_Variable"/> to duplicate. </param>
        public Gh_Variable(Gh_Variable gh_Variable)
        {
            Value = gh_Variable.Value;
        }

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Variable" /> class from a <see cref="GP.Variable"/>.
        /// </summary>
        /// <param name="variable"> <see cref="GP.Variable"/> for the energies and constraints. </param>
        public Gh_Variable(GP.Variable variable)
        {
            Value = variable;
        }

        #endregion


        #region Override : GH_Goo<>

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.IsValid"/>
        public override bool IsValid
        {
            get { return !(Value is null); }
        }

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.TypeDescription"/>
        public override string TypeDescription { get { return string.Format($"Grasshopper type containing a {typeof(GP.Variable)}."); } }

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.TypeName"/>
        public override string TypeName { get { return nameof(Gh_Variable); } }


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.Duplicate"/>
        public override GH_Types.IGH_Goo Duplicate()
        {
            return new Gh_Variable(this);
        }


        /// <inheritdoc cref="GH_Types.GH_Goo{T}.CastFrom(object)"/>
        public override bool CastFrom(object source)
        {
            if (source == null) { return false; }

            var type = source.GetType();

            // ----- BRIDGES Objects ----- //

            // Cast a GP.Variable to a Gh_Variable
            if (typeof(GP.Variable).IsAssignableFrom(type))
            {
                Value = (GP.Variable)source;

                return true;
            }

            // ----- Otherwise ----- //

            return false;
        }

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.CastTo{Q}(ref Q)"/>
        public override bool CastTo<T>(ref T target)
        {

            // ----- BRIDGES Objects ----- //

            // Casts a Gh_Variable to a GP.Variable
            if (typeof(T).IsAssignableFrom(typeof(GP.Variable)))
            {
                object variable = Value;
                target = (T)variable;

                return true;
            }

            // ----- Otherwise ----- //

            return false;
        }

        #endregion

        #region Override : Object

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.ToString"/>
        public override string ToString()
        {
            string text = "Variable: [" + Value[0];
            for (int i = 1; i < Value.Dimension; i++)
            {
                text += ", " + Value[i];
            }
            text += "]";

            return text;
        }

        #endregion
    }
}
