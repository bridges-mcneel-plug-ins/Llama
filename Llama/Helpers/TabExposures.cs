using System;

using GH_Kernel = Grasshopper.Kernel;


namespace Llama.Helpers
{
    /// <summary>
    /// Exposure of the current subcategory tabs.
    /// </summary>
    enum TabExposure
    {
        Mesh = GH_Kernel.GH_Exposure.primary,
        Display = GH_Kernel.GH_Exposure.secondary
    }
}
 