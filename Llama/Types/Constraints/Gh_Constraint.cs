using System;

using GP = BRIDGES.Solvers.GuidedProjection;

using GH_Types = Grasshopper.Kernel.Types;


namespace Llama.Types.Constraints
{
    /// <summary>
    /// Class defining a grasshopper type, which represents an energy.
    /// </summary>
    public class Gh_Constraint : GH_Types.GH_Goo<GP.Constraint>
    {
        #region Constructors

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Constraint" /> class.
        /// </summary>
        public Gh_Constraint() { /* Do Nothing */ }

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Constraint" /> class from another <see cref="Gh_Constraint"/>.
        /// </summary>
        /// <param name="gh_Constraint"> <see cref="Gh_Constraint"/> to duplicate. </param>
        public Gh_Constraint(Gh_Constraint gh_Constraint)
        {
            Value = gh_Constraint.Value;
        }

        /// <summary>
        /// Initialises a new instance of <see cref= "Gh_Constraint" /> class from another <see cref="Gh_Constraint"/>.
        /// </summary>
        /// <param name="constraint"> <see cref="GP.Constraint"/> to store. </param>
        public Gh_Constraint(GP.Constraint constraint)
        {
            Value = constraint;
        }

        #endregion


        #region Override : GH_Goo<>

        // ---------- Properties ---------- //

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.IsValid"/>
        public override bool IsValid => !(Value is null);

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.TypeDescription"/>
        public override string TypeDescription { get { return string.Format($"Grasshopper type containing a {typeof(GP.Constraint)}."); } }

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.TypeName"/>
        public override string TypeName { get { return nameof(Gh_Constraint); } }


        // ---------- Methods ---------- //

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.Duplicate"/>
        public override GH_Types.IGH_Goo Duplicate()
        {
            return new Gh_Constraint(this);
        }


        /// <inheritdoc cref="GH_Types.GH_Goo{T}.CastFrom(object)"/>
        public override bool CastFrom(object source)
        {
            if (source == null) { return false; }

            var type = source.GetType();

            // ----- BRIDGES Objects ----- //

            // Cast a GP.Constraint to a Gh_Constraint
            if (typeof(GP.Constraint).IsAssignableFrom(type))
            {
                Value = (GP.Constraint)source;

                return true;
            }


            // ----- Otherwise ----- //

            return false;
        }

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.CastTo{Q}(ref Q)"/>
        public override bool CastTo<T>(ref T target)
        {
            // ----- BRIDGES Objects ----- //

            // Casts a Gh_Constraint to a GP.Constraint
            if (typeof(T).IsAssignableFrom(typeof(GP.Constraint)))
            {
                object constraint = Value;
                target = (T)constraint;

                return true;
            }

            // ----- Otherwise ----- //

            return false;
        }

        #endregion

        #region Override : Object

        /// <inheritdoc cref="GH_Types.GH_Goo{T}.ToString"/>
        public override string ToString() => $"Constraint (T:{Value.Type})";

        #endregion
    }

    /*
        /// <summary>
        /// Class defining a grasshopper type, which represents an energy.
        /// </summary>
        public class Gh_Constraint : GH_Types.IGH_Goo
        {
            #region Properties

            /// <inheritdoc cref="GH_Types.IGH_Goo.IsValid"/>
            public bool IsValid => true;

            /// <inheritdoc cref="GH_Types.IGH_Goo.IsValidWhyNot"/>
            public string IsValidWhyNot
            {
                get
                {
                    if (IsValid) { return string.Empty; }

                    return $"This {TypeName} is not valid, but I don't know why.";
                }
            }

            /// <inheritdoc cref="GH_Types.IGH_Goo.TypeName"/>
            public string TypeName => nameof(Gh_Constraint);

            /// <inheritdoc cref="GH_Types.IGH_Goo.TypeDescription"/>
            public string TypeDescription => string.Format($"Grasshopper type representing a quadratic constraint for the Guided Projection Algorithm.");


            /// <summary>
            /// Type of the quadratic constraint.
            /// </summary>
            public GP.Interfaces.IQuadraticConstraintType Type { get; private set; }

            /// <summary>
            /// Variables for the quadratic constraint.
            /// </summary>
            public List<Gh_Variable> Variables { get; private set; }

            /// <summary>
            /// Weight of the quadratic constraint.
            /// </summary>
            public double Weight { get; private set; }

            #endregion

            #region Constructors

            /// <summary>
            /// Initialises a new instance of <see cref= "Gh_Constraint" /> class from another <see cref="Gh_Constraint"/>.
            /// </summary>
            /// <param name="gh_QuadraticConstraint"> <see cref="Gh_Constraint"/> to duplicate. </param>
            public Gh_Constraint(Gh_Constraint gh_QuadraticConstraint)
            {
                this.Type = gh_QuadraticConstraint.Type;
                this.Variables = new List<Gh_Variable>(gh_QuadraticConstraint.Variables);

                this.Weight = gh_QuadraticConstraint.Weight;
            }

            /// <summary>
            /// Initialises a new instance of <see cref= "Gh_Constraint" /> class from its name, id, components and variable dimension.
            /// </summary>
            /// <param name="type"> Type of the quadratic constraint. </param>
            /// <param name="variables"> Variables for the quadratic constraint. </param>
            /// <param name="weight"> Weight of the quadratic constraint. </param>
            public Gh_Constraint(GP.Interfaces.IQuadraticConstraintType type, List<Gh_Variable> variables, double weight)
            {
                this.Type = type;
                this.Variables = new List<Gh_Variable>(variables);

                this.Weight = weight;
            }

            #endregion

            #region Public Methods

            /// <inheritdoc cref="GH_Types.IGH_Goo.Duplicate()"/>
            public GH_Types.IGH_Goo Duplicate() => new Gh_Constraint(this);


            /// <inheritdoc cref="GH_Types.IGH_Goo.ScriptVariable()"/>
            public object ScriptVariable() => false;


            /// <inheritdoc cref="GH_Types.IGH_Goo.CastFrom(object)"/>
            public bool CastFrom(object source)
            {
                return false;
            }

            /// <inheritdoc cref="GH_Types.IGH_Goo.CastTo{T}(out T)"/>
            public bool CastTo<T>(out T target)
            {
                target = default;

                return false;
            }


            /// <inheritdoc cref="GH_Types.IGH_Goo.EmitProxy()"/>
            public GH_Types.IGH_GooProxy EmitProxy() => null;


            /// <inheritdoc cref="GH_IO.GH_ISerializable.Read(GH_IO.Serialization.GH_IReader)"/>
            public bool Read(GH_IO.Serialization.GH_IReader reader) => true;

            /// <inheritdoc cref="GH_IO.GH_ISerializable.Write(GH_IO.Serialization.GH_IWriter)"/>
            public bool Write(GH_IO.Serialization.GH_IWriter writer) => true;

            #endregion


            #region Override : Object

            /// <inheritdoc cref="GH_Types.IGH_Goo.ToString()"/>
            public override string ToString() => Type.ToString();

            #endregion
        }
    */

}
