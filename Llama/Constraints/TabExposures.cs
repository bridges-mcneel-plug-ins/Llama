using System;

using GH_Kernel = Grasshopper.Kernel;


namespace Llama.Constraints
{
    /// <summary>
    /// Exposure of the current subcategory tabs.
    /// </summary>
    enum TabExposure
    {
        Numeric = GH_Kernel.GH_Exposure.primary,
        Vector = GH_Kernel.GH_Exposure.secondary,
        Segment = GH_Kernel.GH_Exposure.tertiary,
        Others = GH_Kernel.GH_Exposure.octonary
    }
}
 