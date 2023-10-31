using System;
using System.Collections.Generic;


namespace Llama
{
    /// <summary>
    /// Name of the sub categories. The value is used to sort them of the grasshopper ribbon.
    /// </summary>
    internal enum SubCategory : byte
    {
        Parameters = 0,
        Variables = 1,
        Energies = 2,
        Constraints = 3,
        Solvers = 4
    }

    /// <summary>
    /// Settings for the grasshopper library.
    /// </summary>
    internal static class Settings
    {
        #region Static Properties

        /// <summary>
        /// Name of this category.
        /// </summary>
        public static string CategoryFullName = "BRIDGES - Llama";

        /// <summary>
        /// Nickname of this category.
        /// </summary>
        public static string CategoryName = "Llama";

        /// <summary>
        /// Name of the subcategories.
        /// </summary>
        /// <remarks> 
        /// Tab spaces are accounted for while sorting the subcategories but are trimmed before display (i.e. the more tab spaces, the earlier it will be). 
        /// Therefore, tab spaces are added to arrange the subcategories on the grasshopper ribbon.
        /// </remarks>
        public static Dictionary<SubCategory, string> SubCategoryName = InitialiseSubCategoryName();

        #endregion

        #region Static Methods

        /// <summary>
        /// Initialises <see cref="SubCategoryName"/> with the subcategory names with an appropriate number of spaces.
        /// </summary>
        /// <remarks> 
        /// Tab spaces are accounted for while sorting the subcategories but are trimmed before display (i.e. the more tab spaces, the earlier it will be). 
        /// Therefore, tab spaces are added to arrange the subcategories on the grasshopper ribbon.
        /// </remarks>
        /// <returns> The dictionary with the subcategory names. </returns>
        private static Dictionary<SubCategory, string> InitialiseSubCategoryName()
        {
            Dictionary<SubCategory, string> dictionary = new Dictionary<SubCategory, string>();

            byte max = 0;
            foreach (SubCategory subcatergory in (SubCategory[])Enum.GetValues(typeof(SubCategory)))
            {
                if (max < (byte)subcatergory) { max = (byte)subcatergory; }
            }
            foreach (SubCategory subcatergory in (SubCategory[])Enum.GetValues(typeof(SubCategory)))
            {
                dictionary.Add(subcatergory, new string('\t', max - (byte)subcatergory) + subcatergory.ToString());
            }

            return dictionary;
        }

        #endregion
    }
}
