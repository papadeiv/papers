"""
Color palettes used in the plots.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 26-09-2025
"""

# Catpuccinn frappe'
CtpMauve = colorant"rgb(202,158,230)"
CtpTeal = colorant"rgb(129, 200, 190)"
CtpGreen = colorant"rgb(166, 209, 137)"
CtpBlue = colorant"rgb(140, 170, 238)"
CtpRed = colorant"rgb(231, 130, 132)"
CtpYellow = colorant"rgb(229,200,144)"
CtpWhite = colorant"rgb(198,208,245)"
CtpGray = colorant"rgb(98,104,128)"

# Export the palette
for name in names(@__MODULE__; all = true)
        @eval export $name
end
