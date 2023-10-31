using System;

using GH_Kernel = Grasshopper.Kernel;


namespace Llama.Models
{
    /// <summary>
    /// Exposure of the current subcategory tabs.
    /// </summary>
    enum TabExposure
    {
        Variable = GH_Kernel.GH_Exposure.primary,
        Energy = GH_Kernel.GH_Exposure.secondary,
        Constraint = GH_Kernel.GH_Exposure.tertiary,
        Solver = GH_Kernel.GH_Exposure.quarternary,
        Nexorade = GH_Kernel.GH_Exposure.quinary,
    }
}
