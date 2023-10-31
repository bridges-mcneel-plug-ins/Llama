using System;

using GH_Kernel = Grasshopper.Kernel;


namespace Llama.Parameters
{
    /// <summary>
    /// Exposure of the current subcategory tabs.
    /// </summary>
    enum TabExposure
    {
        Variables = GH_Kernel.GH_Exposure.primary,
        Energies = GH_Kernel.GH_Exposure.secondary,
        Constraints = GH_Kernel.GH_Exposure.tertiary,
        Models = GH_Kernel.GH_Exposure.quarternary
    }
}
