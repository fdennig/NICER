#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plots
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set User
user = "joshua" # or "francis" as the case may be

# Load Data
regime_select = 2
  # 0 = no randomization (just uses means)
  # 1 = total randomization

  # 2 = High initial TFP growth vs. Low initial TFP growth
  # 3 = High initial decarbonization rate vs. Low initial decarbonization rate
  # 4 = High elasticity of income wrt damage vs. Low elasticity of income wrt damage
  # 5 = High climate sensitivity vs. Low climate sensitivity
  # 6 = High atmosphere to upper ocean transfer coefficient vs. Low atmosphere to upper ocean transfer coefficient
  # 7 = High initial world backstop price vs. Low initial world backstop price
  # 8 = High T7 coefficient vs. Low T7 coefficient

  # 9 = Deciles of TFP (all else fixed at means) (nsample = 10)
  # 10 = Deciles of decarbonization rates (all else fixed at means) (nsample = 10)
  # 11 = Deciles - High TFP and High decarb vs. Low TFP and Low decarb (technology spillovers?)
  # 12 = Deciles - High TFP and Low decarb vs. Low TFP and High decarb (substitutable tech?)
  # 13 = Deciles of elasticity of income wrt damage (ee)
  # 14 = Deciles - High TFP and High ee vs. Low TFP and Low ee
  # 15 = Deciles - High TFP and Low ee vs. Low TFP and High ee

# Set directory based on User
if user == "francis"
  folder = "/Users/francis/Dropbox/ARBEIT/aRESEARCH/"
elseif user == "joshua"
  folder = "/Users/joshuabernstein/Dropbox"
elseif user == "marc"
  folder = "\Users\mfleur\Dropbox\RICE_for_Simon (1)\Julia"
else error("wrong user")
end

dt = string(regime_select)
using JLD
Data = load("$folder/NICE Julia/Outputs/SS_$dt.jld")
# To call data from within the Data object, the syntax is Data["SS"][i].'variable name' where i refers to the choice of lm (1->0, 2->4, 3->tm), and the 'variable name' is one of the following:
#   regime
#   nsample
#   Tm
#   tm
#   lm
#   Regions
#   taxes_1
#   taxes_2
#   EWelfare
#   c
#   K
#   T
#   E
#   M
#   mu
#   lam
#   D
#   AD
#   Y
#   Q
#   rho
#   eta
#   nu
#   PP

# create the x axis "time"
tm = Data["SS"][1].tm
time_tm = (10.*[1:tm]).+2015 # years for taxes

# Get y-axis data from Data object
# e.g. get taxes from Data object
Y_Array = [Data["SS"][1].taxes_1 Data["SS"][1].taxes_2 Data["SS"][2].taxes_1 Data["SS"][2].taxes_2 Data["SS"][3].taxes_1]

# Create DataFrame for plots
Cats = [repmat(["lm = 0, Higher"],tm),repmat(["lm = 0, Lower"],tm),repmat(["lm = 4, Higher"],tm),repmat(["lm = 4, Lower"],tm),repmat(["lm = tm"],tm)]
using DataFrames
taxes = DataFrame(Year = repmat(time_tm,5), Tax = vec(Y_Array), Category = Cats)

# Plot the tax paths over time
# fnames = ["TFP" "Decarb" "Elas_inc" "Clim_Sens" "Atmos_Oce" "Backstop" "T7_Coeff" "Deciles TFP" "Deciles decarb" "Deciles High TFP and decarb" "Deciles High TFP low decarb" "Deciles ee" "Deciles High TFP and ee" "Deciles High TFP Low ee"]
filenm = string(regime_select) # filename
using Gadfly, Cairo
Tax_Plot = plot(taxes,x="Year",y="Tax",color="Category",Scale.x_continuous(minvalue=time_tm[1], maxvalue=time_tm[tm]),Guide.title(""),Geom.line,Guide.colorkey("")) # Scale.color_discrete_manual("dark green","light green", "dark blue","light blue","magenta") ,"blue","green","yellow","black"
  #draw(PDF("$folder/NICE Julia/Plots/5_line_$filenm.pdf", 6inch, 4.5inch),Tax_Plot)

