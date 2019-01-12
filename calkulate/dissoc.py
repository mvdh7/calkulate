from numpy import exp, log


# ===== COMBINED KX LISTS =====================================================

# --- DSC07 best practice, Free scale -----------------------------------------

def KX_F(Tk,S,ST,FT):
    
    # Evaluate dissociation coefficients
    KC1_T,KC2_T       = KC_T_LDK00  (Tk,S) # Total
    KB_T              = KB_T_D90a   (Tk,S) # Total
    KH2O_T            = KH2O_T_DSC07(Tk,S) # Total
    KHSO4_F           = KHSO4_F_D90b(Tk,S) # Free
    KHF_T             = KHF_T_PF87  (Tk,S) # Total
    KP1_T,KP2_T,KP3_T = KP_T_DSC07  (Tk,S) # Total
    KSi_T             = KSi_T_M95   (Tk,S) # Total
    
    # Get pH scale conversion factors
    T2F = 1 / (1 + ST/KHSO4_F)
    
    # Convert everything to Free pH scale
    KC1_F  = KC1_T  * T2F
    KC2_F  = KC2_T  * T2F
    KB_F   = KB_T   * T2F
    KH2O_F = KH2O_T * T2F
    KHF_F  = KHF_T  * T2F
    KP1_F  = KP1_T  * T2F
    KP2_F  = KP2_T  * T2F
    KP3_F  = KP3_T  * T2F
    KSi_F  = KSi_T  * T2F
    
    return [KC1_F,KC2_F, KB_F, KH2O_F, KHSO4_F, KHF_F, KP1_F,KP2_F,KP3_F,
            KSi_F]


# ===== IONIC STRENGTH ========================================================


def Istr(S):
    return 19.924 * S / (1000 - 1.005 * S)


# ===== CARBONIC ACID =========================================================

# --- Lueker et al. (2000) ----------------------------------------------------
#
# Mar Chem 70(1-3), 105-119, doi:10.1016/S0304-4203(00)00022-0
#
# Total pH scale
#
#  2 < T < 35 degC
# 19 < S < 43

def KC_T_LDK00(Tk,S):

    # LDK00 Eq. (16)
    pK_T_C1 = 3633.86      /     Tk  \
            -   61.2172              \
            +    9.6777    * log(Tk) \
            -    0.011555  * S       \
            +    0.0001152 * S**2

    # LDK00 Eq. (17)
    pK_T_C2 =  471.78      /     Tk  \
            +   25.929               \
            -    3.16967   * log(Tk) \
            -    0.01781   * S       \
            +    0.0001122 * S**2

    return 10**-pK_T_C1, 10**-pK_T_C2


# ===== BORIC ACID ============================================================

# --- Dickson (1990a) ---------------------------------------------------------
#
# Deep-Sea Res Pt A 37(5), 755-766, doi:10.1016/0198-0149(90)90004-F
#
# Total pH scale
# 
#  0 < T < 45 degC
#  5 < S < 45
    
def KB_T_D90a(Tk,S):

    # D90a Eq. (23)
    ln_KB =   ( - 8966.90                        \
                - 2890.53   * S**0.5             \
                -   77.942  * S                  \
                +    1.728  * S**1.5             \
                -    0.0996 * S**2   ) /     Tk  \
            + 148.0248                           \
            + 137.1942  * S**0.5                 \
            +   1.62142 * S                      \
            - (     24.4344                      \
                +   25.085  * S**0.5             \
                +    0.2474 * S      ) * log(Tk) \
            +   0.053105    * S**0.5   *     Tk

    return exp(ln_KB)


# ===== WATER =================================================================

# --- Dickson et al. (2007) ---------------------------------------------------
#
# PICES Special Publication 3
#
# Total pH scale

def KH2O_T_DSC07(Tk,S):

    ln_KH2O =     148.9652                        \
              - 13847.26    /     Tk              \
              -    23.6521  * log(Tk)             \
              + (   118.67   /     Tk             \
                  -   5.977                       \
                  +   1.0495 * log(Tk) ) * S**0.5 \
              - 0.01615                 * S

    return exp(ln_KH2O)


# ===== BISULFATE =============================================================

# --- Dickson (1990b) ---------------------------------------------------------
#
# J Chem Thermodyn 22(2), 113-127, doi:10.1016/0021-9614(90)90074-Z
#
# Free pH scale
# 
#  0 < T < 45 degC
#  5 < S < 45

def KHSO4_F_D90b(Tk,S):

    # Ionic strength
    I = Istr(S)

    # D90b Eqs. (22) & (23)
    ln_KHSO4 = - 4276.1   /     Tk                  \
               +  141.328                           \
               -   23.093 * log(Tk)                 \
               + ( - 13856     /     Tk             \
                   +   324.57                       \
                   -    47.986 * log(Tk) ) * I**0.5 \
               + (   35474     /     Tk             \
                   -   771.54                       \
                   +   114.723 * log(Tk) ) * I      \
               - (    2698     /     Tk)   * I**1.5 \
               + (    1776     /     Tk)   * I**2   \
               + log(1 - 0.001005 * S)

    return exp(ln_KHSO4)


# ===== HYDROGEN FLUORIDE =====================================================

# --- Perez & Fraga (1987) ----------------------------------------------------
#
# Mar Chem 21(2), 161-168, doi:10.1016/0304-4203(87)90036-3
#
# Total pH scale
# 
#  9 < T < 33 degC
# 10 < S < 40

def KHF_T_PF87(Tk,S):
    
    ln_KHF = - ( - 874     / Tk        \
                 -   0.111 * S **0.5   \
                 +   9.68            )
             
    return exp(ln_KHF)
             

# --- Dickson & Riley (1979) --------------------------------------------------
#
# Mar Chem 7(2), 101-109, doi:10.1016/0304-4203(79)90002-1
#
# Free pH scale
# 
#  5 < T < 35 degC
# 10 < S < 48

def KHF_F_DR79(Tk,S):
    
    # Ionic strength
    I = Istr(S)
    
    # Evaluate HF dissociation constant
    ln_KF =   1590.2   / Tk       \
            -   12.641            \
            +    1.525 * I **0.5  \
            + log(1 - 0.001005*S)
            
    return exp(ln_KF)
    

# ===== PHOSPHORIC ACID =======================================================

# --- Dickson et al. (2007) ---------------------------------------------------
#
# PICES Special Publication 3
#
# Total pH scale

def KP_T_DSC07(Tk,S):
    
    # KP1 = [H+][H2PO4-]/[H3PO4]
    ln_KP1 = - 4576.752 /     Tk            \
            +  115.525                     \
            -   18.453 * log(Tk)           \
            + (- 106.736    / Tk           \
               +   0.69171      ) * S**0.5 \
            + (-   0.65643 / Tk            \
               -   0.01844      ) * S

    # KP2 = [H+][HPO42-]/[H2PO4-]
    ln_KP2 = - 8814.715 /     Tk           \
            +  172.0883                   \
            -   27.927 * log(Tk)          \
            + (- 160.34    / Tk           \
               +   1.3566      ) * S**0.5 \
            + (    0.37335 / Tk           \
               -   0.05778     ) * S

    # KP3 = [H+][PO43-]/[HPO42-]
    ln_KP3 = - 3070.75 / Tk               \
            - 18.141                     \
            + (  17.27039 / Tk           \
               +  2.81197     ) * S**0.5 \
            + (- 44.99486 / Tk           \
               -  0.09984     ) * S

    return exp(ln_KP1), exp(ln_KP2), exp(ln_KP3)


# ===== SILICIC ACID ==========================================================
    
# --- Millero (1995) ----------------------------------------------------------
#
#
#
# via Dickson et al. (2007)

def KSi_T_M95(Tk,S):
    
    # Ionic strength
    I = 19.924 * S / (1000 - 1.005 * S)
    
    #  KSi = [SiO((OH)3)-] [H+] / [Si(OH)4]
    ln_KSi = - 8904.2   /     Tk   \
             +  117.385            \
             -   19.334 * log(Tk)  \
             + ( - 458.79    / Tk   \
                 +   3.5913       ) * I**0.5 \
             + (   188.74    / Tk            \
                 -   1.5998       ) * I      \
             + ( -  12.1652  / Tk            \
                 +   0.07871      ) * I**2   \
             + log(1 - 0.001005 * S)
             
    return exp(ln_KSi)


# ===== AMMONIUM ==============================================================
    
# --- Bell et al. (2008) ------------------------------------------------------
#    
# Environ Chem 5, 258, doi:10.1071/EN07032_CO
#
# pH scale unclear

def KNH4_X_BJJL08(Tk,S):
    
    # BJJL08 Eq. (3)
    pKNH4 = 10.0423         \
          -  0.0315536 * Tk \
          +  0.003071  * S
          
    return 10**-pKNH4


# ===== 2-AMINOPYRIDINE =======================================================

# --- Bates & Erickson (1986) -------------------------------------------------
#
# J Solution Chem 15(11), 891-901, doi:10.1007/BF00646030
#
# Stoichiometric dissociation constant
#
# Seawater pH scale
# 
#  5 < T < 40 degC
# 30 < S < 40

def K2AMP_S_BE86(Tk,S):

    # BE86 Eq. (10)
    pK_S_2AMP =   2498.31   /     Tk        \
                -   15.3274                 \
                +    2.4050 * log(Tk)       \
                + (   0.012929              \
                    - 2.9417e-5 * Tk  ) * S

    return 10**-pK_S_2AMP
