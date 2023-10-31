using System;

using GH_Kernel = Grasshopper.Kernel;


namespace Llama.Energies
{
    /// <summary>
    /// Exposure of the current subcategory tabs.
    /// </summary>
    enum TabExposure
    {
        Numeric = GH_Kernel.GH_Exposure.primary,
        Vector = GH_Kernel.GH_Exposure.secondary,
        Segment = GH_Kernel.GH_Exposure.tertiary,
        Face = GH_Kernel.GH_Exposure.quarternary,
        Others = GH_Kernel.GH_Exposure.octonary,
    }
}
 