--------------------------------------------------------------------------------
-- Quanty input file generated using Crispy. If you use this file please cite
-- the following reference: 10.5281/zenodo.1008184.
--
-- elements: 3d transition metals
-- symmetry: Oh
-- experiment: XAS
-- edge: L2,3 (2p)
--------------------------------------------------------------------------------
Verbosity(0x0000)

--------------------------------------------------------------------------------
-- Initialize the Hamiltonians.
--------------------------------------------------------------------------------
H_i = 0
H_f = 0

--------------------------------------------------------------------------------
-- Toggle the Hamiltonian terms.
--------------------------------------------------------------------------------
H_atomic              = 1
H_cf                  = 1
H_3d_Ld_hybridization = 1
H_magnetic_field      = 1
H_exchange_field      = 0

--------------------------------------------------------------------------------
-- Define the number of electrons, shells, etc.
--------------------------------------------------------------------------------
NBosons = 0
NFermions = 16

NElectrons_2p = 6
NElectrons_3d = 8

IndexDn_2p = {0, 2, 4}
IndexUp_2p = {1, 3, 5}
IndexDn_3d = {6, 8, 10, 12, 14}
IndexUp_3d = {7, 9, 11, 13, 15}

if H_3d_Ld_hybridization == 1 then
    NFermions = 26

    NElectrons_Ld = 10

    IndexDn_Ld = {16, 18, 20, 22, 24}
    IndexUp_Ld = {17, 19, 21, 23, 25}
end

--------------------------------------------------------------------------------
-- Define the atomic term.
--------------------------------------------------------------------------------
N_2p = NewOperator('Number', NFermions, IndexUp_2p, IndexUp_2p, {1, 1, 1})
     + NewOperator('Number', NFermions, IndexDn_2p, IndexDn_2p, {1, 1, 1})

N_3d = NewOperator('Number', NFermions, IndexUp_3d, IndexUp_3d, {1, 1, 1, 1, 1})
     + NewOperator('Number', NFermions, IndexDn_3d, IndexDn_3d, {1, 1, 1, 1, 1})

if H_atomic == 1 then
    F0_3d_3d = NewOperator('U', NFermions, IndexUp_3d, IndexDn_3d, {1, 0, 0})
    F2_3d_3d = NewOperator('U', NFermions, IndexUp_3d, IndexDn_3d, {0, 1, 0})
    F4_3d_3d = NewOperator('U', NFermions, IndexUp_3d, IndexDn_3d, {0, 0, 1})

    F0_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {1, 0}, {0, 0})
    F2_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 1}, {0, 0})
    G1_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 0}, {1, 0})
    G3_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 0}, {0, 1})

    U_3d_3d_i  = 7.3 * 1.0
    F2_3d_3d_i = 12.234 * 0.8
    F4_3d_3d_i = 7.598 * 0.8
    F0_3d_3d_i = U_3d_3d_i + 2 / 63 * F2_3d_3d_i + 2 / 63 * F4_3d_3d_i

    U_3d_3d_f  = 7.3 * 1.0
    F2_3d_3d_f = 13.0067 * 0.8
    F4_3d_3d_f = 8.0845 * 0.8
    F0_3d_3d_f = U_3d_3d_f + 2 / 63 * F2_3d_3d_f + 2 / 63 * F4_3d_3d_f
    U_2p_3d_f  = 8.3 * 1.0
    F2_2p_3d_f = 7.721 * 0.8
    G1_2p_3d_f = 5.787 * 0.8
    G3_2p_3d_f = 3.291 * 0.8
    F0_2p_3d_f = U_2p_3d_f + 1 / 15 * G1_2p_3d_f + 3 / 70 * G3_2p_3d_f

    H_i = H_i
        + F0_3d_3d_i * F0_3d_3d
        + F2_3d_3d_i * F2_3d_3d
        + F4_3d_3d_i * F4_3d_3d

    H_f = H_f
        + F0_3d_3d_f * F0_3d_3d
        + F2_3d_3d_f * F2_3d_3d
        + F4_3d_3d_f * F4_3d_3d
        + F0_2p_3d_f * F0_2p_3d
        + F2_2p_3d_f * F2_2p_3d
        + G1_2p_3d_f * G1_2p_3d
        + G3_2p_3d_f * G3_2p_3d

    ldots_3d = NewOperator('ldots', NFermions, IndexUp_3d, IndexDn_3d)

    ldots_2p = NewOperator('ldots', NFermions, IndexUp_2p, IndexDn_2p)

    zeta_3d_i = 0.083 * 1.0

    zeta_3d_f = 0.102 * 1.0
    zeta_2p_f = 11.507 * 1.0

    H_i = H_i
        + zeta_3d_i * ldots_3d

    H_f = H_f
        + zeta_3d_f * ldots_3d
        + zeta_2p_f * ldots_2p
end

--------------------------------------------------------------------------------
-- Define the crystal field term.
--------------------------------------------------------------------------------
if H_cf == 1 then
    -- PotentialExpandedOnClm('D4h', 2, {Ea1g, Eb1g, Eb2g, Eeg})
    Dq_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, { 6,  6, -4, -4}))
    Ds_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, {-2,  2,  2, -1}))
    Dt_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, {-6, -1, -1,  4}))

    Dq_3d_i = 0.15
    Ds_3d_i = 0.02
    Dt_3d_i = 0.01

    Dq_3d_f = 0.15
    Ds_3d_f = 0.02
    Dt_3d_f = 0.01

    H_i = H_i
        + Dq_3d_i * Dq_3d
        + Ds_3d_i * Ds_3d
        + Dt_3d_i * Dt_3d

    H_f = H_f
        + Dq_3d_f * Dq_3d
        + Ds_3d_f * Ds_3d
        + Dt_3d_f * Dt_3d
end

--------------------------------------------------------------------------------
-- Define the 3d-Ld hybridization term.
--------------------------------------------------------------------------------
if H_3d_Ld_hybridization == 1 then
    N_Ld = NewOperator('Number', NFermions, IndexUp_Ld, IndexUp_Ld, {1, 1, 1, 1, 1})
         + NewOperator('Number', NFermions, IndexDn_Ld, IndexDn_Ld, {1, 1, 1, 1, 1})

    Delta_3d_Ld_i = 5.3
    e_3d_i  = (10 * Delta_3d_Ld_i - NElectrons_3d * (19 + NElectrons_3d) * U_3d_3d_i / 2) / (10 + NElectrons_3d)
    e_Ld_i  = NElectrons_3d * ((1 + NElectrons_3d) * U_3d_3d_i / 2 - Delta_3d_Ld_i) / (10 + NElectrons_3d)

    Delta_3d_Ld_f = 5.3
    e_3d_f = (10 * Delta_3d_Ld_f - NElectrons_3d * (31 + NElectrons_3d) * U_3d_3d_f / 2 - 90 * U_2p_3d_f) / (16 + NElectrons_3d)
    e_2p_f = (10 * Delta_3d_Ld_f + (1 + NElectrons_3d) * (NElectrons_3d * U_3d_3d_f / 2 - (10 + NElectrons_3d) * U_2p_3d_f)) / (16 + NElectrons_3d)
    e_Ld_f = ((1 + NElectrons_3d) * (NElectrons_3d * U_3d_3d_f / 2 + 6 * U_2p_3d_f) - (6 + NElectrons_3d) * Delta_3d_Ld_f) / (16 + NElectrons_3d)

    H_i = H_i
        + e_3d_i * N_3d
        + e_Ld_i * N_Ld

    H_f = H_f
        + e_3d_f * N_3d
        + e_2p_f * N_2p
        + e_Ld_f * N_Ld

    Dq_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, { 6,  6, -4, -4}))
    Ds_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, {-2,  2,  2, -1}))
    Dt_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, {-6, -1, -1,  4}))

    Va1g_3d_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, {1, 0, 0, 0}))
               + NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, {1, 0, 0, 0}))

    Vb1g_3d_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, {0, 1, 0, 0}))
               + NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, {0, 1, 0, 0}))

    Vb2g_3d_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, {0, 0, 1, 0}))
               + NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, {0, 0, 1, 0}))

    Veg_3d_Ld = NewOperator('CF', NFermions, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm('D4h', 2, {0, 0, 0, 1}))
              + NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, PotentialExpandedOnClm('D4h', 2, {0, 0, 0, 1}))

    Dq_Ld_i = 0.0
    Ds_Ld_i = 0.0
    Dt_Ld_i = 0.0
    Va1g_3d_Ld_i = 3.3
    Vb1g_3d_Ld_i = 4.4
    Vb2g_3d_Ld_i = 2.2
    Veg_3d_Ld_i  = 1.1

    Dq_Ld_f = 0.0
    Ds_Ld_f = 0.0
    Dt_Ld_f = 0.0
    Va1g_3d_Ld_f = 3.3
    Vb1g_3d_Ld_f = 4.4
    Vb2g_3d_Ld_f = 2.2
    Veg_3d_Ld_f  = 1.1

    H_i = H_i
        + Dq_Ld_i      * Dq_Ld
        + Ds_Ld_i      * Ds_Ld
        + Dt_Ld_i      * Dt_Ld
        + Va1g_3d_Ld_i * Va1g_3d_Ld
        + Vb1g_3d_Ld_i * Vb1g_3d_Ld
        + Vb2g_3d_Ld_i * Vb2g_3d_Ld
        + Veg_3d_Ld_i  * Veg_3d_Ld

    H_f = H_f
        + Dq_Ld_f      * Dq_Ld
        + Ds_Ld_f      * Ds_Ld
        + Dt_Ld_f      * Dt_Ld
        + Va1g_3d_Ld_f * Va1g_3d_Ld
        + Vb1g_3d_Ld_f * Vb1g_3d_Ld
        + Vb2g_3d_Ld_f * Vb2g_3d_Ld
        + Veg_3d_Ld_f  * Veg_3d_Ld
end

--------------------------------------------------------------------------------
-- Define the magnetic field and exchange field terms.
--------------------------------------------------------------------------------
Sx_3d    = NewOperator('Sx'   , NFermions, IndexUp_3d, IndexDn_3d)
Sy_3d    = NewOperator('Sy'   , NFermions, IndexUp_3d, IndexDn_3d)
Sz_3d    = NewOperator('Sz'   , NFermions, IndexUp_3d, IndexDn_3d)
Ssqr_3d  = NewOperator('Ssqr' , NFermions, IndexUp_3d, IndexDn_3d)
Splus_3d = NewOperator('Splus', NFermions, IndexUp_3d, IndexDn_3d)
Smin_3d  = NewOperator('Smin' , NFermions, IndexUp_3d, IndexDn_3d)

Lx_3d    = NewOperator('Lx'   , NFermions, IndexUp_3d, IndexDn_3d)
Ly_3d    = NewOperator('Ly'   , NFermions, IndexUp_3d, IndexDn_3d)
Lz_3d    = NewOperator('Lz'   , NFermions, IndexUp_3d, IndexDn_3d)
Lsqr_3d  = NewOperator('Lsqr' , NFermions, IndexUp_3d, IndexDn_3d)
Lplus_3d = NewOperator('Lplus', NFermions, IndexUp_3d, IndexDn_3d)
Lmin_3d  = NewOperator('Lmin' , NFermions, IndexUp_3d, IndexDn_3d)

Jx_3d    = NewOperator('Jx'   , NFermions, IndexUp_3d, IndexDn_3d)
Jy_3d    = NewOperator('Jy'   , NFermions, IndexUp_3d, IndexDn_3d)
Jz_3d    = NewOperator('Jz'   , NFermions, IndexUp_3d, IndexDn_3d)
Jsqr_3d  = NewOperator('Jsqr' , NFermions, IndexUp_3d, IndexDn_3d)
Jplus_3d = NewOperator('Jplus', NFermions, IndexUp_3d, IndexDn_3d)
Jmin_3d  = NewOperator('Jmin' , NFermions, IndexUp_3d, IndexDn_3d)

Sx = Sx_3d
Sy = Sy_3d
Sz = Sz_3d

Lx = Lx_3d
Ly = Ly_3d
Lz = Lz_3d

Jx = Jx_3d
Jy = Jy_3d
Jz = Jz_3d

Ssqr = Sx * Sx + Sy * Sy + Sz * Sz
Lsqr = Lx * Lx + Ly * Ly + Lz * Lz
Jsqr = Jx * Jx + Jy * Jy + Jz * Jz

if H_magnetic_field == 1 then
    Bx_i = 0.0 -- 0.0 * EnergyUnits.Tesla.value
    By_i = 0.0 -- 0.0 * EnergyUnits.Tesla.value
    Bz_i = 0.001 -- 10.0 * EnergyUnits.Tesla.value

    Bx_f = 0.0 -- 0.0 * EnergyUnits.Tesla.value
    By_f = 0.0 -- 0.0 * EnergyUnits.Tesla.value
    Bz_f = 0.001 -- 10.0 * EnergyUnits.Tesla.value

    H_i = H_i
        + Bx_i * (2 * Sx + Lx)
        + By_i * (2 * Sy + Ly)
        + Bz_i * (2 * Sz + Lz)

    H_f = H_f
        + Bx_f * (2 * Sx + Lx)
        + By_f * (2 * Sy + Ly)
        + Bz_f * (2 * Sz + Lz)
end

if H_exchange_field == 1 then
    Hx_i = 0.0
    Hy_i = 0.0
    Hz_i = 0.002

    Hx_f = 0.0
    Hy_f = 0.0
    Hz_f = 0.002

    H_i = H_i
        + Hx_i * Sx
        + Hy_i * Sy
        + Hz_i * Sz

    H_f = H_f
        + Hx_f * Sx
        + Hy_f * Sy
        + Hz_f * Sz
end

NConfigurations = 2

--------------------------------------------------------------------------------
-- Define the restrictions and set the number of initial states.
--------------------------------------------------------------------------------
InitialRestrictions = {NFermions, NBosons, {'111111 0000000000', NElectrons_2p, NElectrons_2p},
                                           {'000000 1111111111', NElectrons_3d, NElectrons_3d}}

FinalRestrictions = {NFermions, NBosons, {'111111 0000000000', NElectrons_2p - 1, NElectrons_2p - 1},
                                         {'000000 1111111111', NElectrons_3d + 1, NElectrons_3d + 1}}

if H_3d_Ld_hybridization == 1 then
    InitialRestrictions = {NFermions, NBosons, {'111111 0000000000 0000000000', NElectrons_2p, NElectrons_2p},
                                               {'000000 1111111111 0000000000', NElectrons_3d, NElectrons_3d},
                                               {'000000 0000000000 1111111111', NElectrons_Ld, NElectrons_Ld}}

    FinalRestrictions = {NFermions, NBosons, {'111111 0000000000 0000000000', NElectrons_2p - 1, NElectrons_2p - 1},
                                             {'000000 1111111111 0000000000', NElectrons_3d + 1, NElectrons_3d + 1},
                                             {'000000 0000000000 1111111111', NElectrons_Ld, NElectrons_Ld}}

    CalculationRestrictions = {NFermions, NBosons, {'000000 0000000000 1111111111', NElectrons_Ld - (NConfigurations - 1), NElectrons_Ld}}
end

Operators = {H_i, Ssqr, Lsqr, Jsqr, Sz, Lz, Jz, N_2p, N_3d}
header = 'Analysis of the initial Hamiltonian:\n'
header = header .. '==============================================================================================\n'
header = header .. '   i       <E>     <S^2>     <L^2>     <J^2>      <Sz>      <Lz>      <Jz>    <N_2p>    <N_3d>\n'
header = header .. '==============================================================================================\n'
footer = '==============================================================================================\n'

if H_3d_Ld_hybridization == 1 then
    Operators = {H_i, Ssqr, Lsqr, Jsqr, Sz, Lz, Jz, N_2p, N_3d, N_Ld}
    header = 'Analysis of the initial Hamiltonian:\n'
    header = header .. '========================================================================================================\n'
    header = header .. '   i       <E>     <S^2>     <L^2>     <J^2>      <Sz>      <Lz>      <Jz>    <N_2p>    <N_3d>    <N_Ld>\n'
    header = header .. '========================================================================================================\n'
    footer = '========================================================================================================\n'
end

-- Define the temperature.
T = 0.0 * EnergyUnits.Kelvin.value

 -- Approximate machine epsilon.
epsilon = 2.22e-16
Z = 0

NPsis = 16
NPsisAuto = 0

if NPsisAuto == 1 and NPsis ~= 1 then
    NPsis = 1
    NPsisIncrement = 8
    NPsisIsConverged = false
    dZ = {}

    while not NPsisIsConverged do
        if CalculationRestrictions == nil then
            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis)
        else
            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis, {{'restrictions', CalculationRestrictions}})
        end

        if not (type(Psis_i) == 'table') then
            Psis_i = {Psis_i}
        end

        E_gs_i = Psis_i[1] * H_i * Psis_i[1]

        for i, Psi in ipairs(Psis_i) do
            E = Psi * H_i * Psi

            if math.abs(E - E_gs_i) < epsilon then
                dZ[i] = 1
            else
                dZ[i] = math.exp(-(E - E_gs_i) / T)
            end

            Z = Z + dZ[i]

            if (dZ[i] / Z) < math.sqrt(epsilon) then
                i = i - 1
                NPsisIsConverged = true
                NPsis = i
                Psis_i = {unpack(Psis_i, 1, i)}
                dZ = {unpack(dZ, 1, i)}
                break
            end
        end

        if NPsisIsConverged then
            break
        else
            NPsis = NPsis + NPsisIncrement
        end
    end
else
        if CalculationRestrictions == nil then
            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis)
        else
            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis, {{'restrictions', CalculationRestrictions}})
        end

    if not (type(Psis_i) == 'table') then
        Psis_i = {Psis_i}
    end
end

io.write(header)
for i, Psi in ipairs(Psis_i) do
    io.write(string.format('%4d', i))
    for j, Operator in ipairs(Operators) do
        io.write(string.format('%10.5f', Complex.Re(Psi * Operator * Psi)))
    end
    io.write('\n')
end
io.write(footer)

--------------------------------------------------------------------------------
-- Define the transition operators.
--------------------------------------------------------------------------------
t = math.sqrt(1/2);

kin = {0, 0, -1}
ein1 = {0, 1, 0}
ein2 = {-1, 0, 0}

Tx_2p_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1, -1, t    }, {1, 1, -t    }})
Ty_2p_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1, -1, t * I}, {1, 1,  t * I}})
Tz_2p_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1,  0, 1    }                })

Tein1_2p_3d = ein1[1] * Tx_2p_3d + ein1[2] * Ty_2p_3d + ein1[3] * Tz_2p_3d
Tein2_2p_3d = ein2[1] * Tx_2p_3d + ein2[2] * Ty_2p_3d + ein2[3] * Tz_2p_3d

Tr_2p_3d =  t * (Tein1_2p_3d - I * Tein2_2p_3d)
Tl_2p_3d = -t * (Tein1_2p_3d + I * Tein2_2p_3d)

--------------------------------------------------------------------------------
-- Calculate and save the spectra.
--------------------------------------------------------------------------------
CalculateIso = 1
CalculateCD  = 1
CalculateLD  = 1

if CalculateIso == 0 and CalculateCD == 0 and CalculateLD == 0 then
    return
end

E_gs_i = Psis_i[1] * H_i * Psis_i[1]

NPsis_f = 16

if CalculationRestrictions == nil then
    Psis_f = Eigensystem(H_f, FinalRestrictions, NPsis_f)
else
    Psis_f = Eigensystem(H_f, FinalRestrictions, NPsis_f, {{'restrictions', CalculationRestrictions}})
end

Operators = {H_f, Ssqr, Lsqr, Jsqr, Sz, Lz, Jz, N_2p, N_3d, N_Ld}
header = 'Analysis of the final Hamiltonian:\n'
header = header .. '========================================================================================================\n'
header = header .. '   i       <E>     <S^2>     <L^2>     <J^2>      <Sz>      <Lz>      <Jz>    <N_2p>    <N_3d>    <N_Ld>\n'
header = header .. '========================================================================================================\n'
footer = '========================================================================================================\n'

io.write('\n')
io.write(header)
for i, Psi in ipairs(Psis_f) do
    io.write(string.format('%4d', i))
    for j, Operator in ipairs(Operators) do
        io.write(string.format('%10.5f', Complex.Re(Psi * Operator * Psi)))
    end
    io.write('\n')
end
io.write(footer)

-- Psis_f = {Psis_f}
E_gs_f = Psis_f[1] * H_f * Psis_f[1]

Eedge1 = 852.7
DeltaE = Eedge1 + E_gs_i - E_gs_f

Emin = 842.7 - DeltaE
Emax = 882.7 - DeltaE
Gamma = 0.1
NE = 3999

Z = 0

Giso = 0

Gr = 0
Gl = 0

Gein1 = 0
Gein2 = 0

io.write(string.format('\nSpectrum calculation for each of the selected states:\n'))
io.write(string.format('===============\n'))
io.write(string.format('   i         dZ\n'))
io.write(string.format('===============\n'))

Psis_i = {Psis_i[1], Psis_i[2], Psis_i[3]}

for i, Psi in ipairs(Psis_i) do
    E = Psi * H_i * Psi

    if math.abs(E - E_gs_i) < epsilon then
        dZ = 1
    else
        dZ = math.exp(-(E - E_gs_i) / T)
    end

    Z = Z + dZ

    io.write(string.format('%4d   %3.2E\n', i, dZ))

    if CalculateIso == 1 then
        for j, Operator in ipairs({Tx_2p_3d, Ty_2p_3d, Tz_2p_3d}) do
            if CalculationRestrictions == nil then
                Giso = Giso + CreateSpectra(H_f, Operator, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}}) * dZ
            else
                Giso = Giso + CreateSpectra(H_f, Operator, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'restrictions', CalculationRestrictions}}) * dZ
            end
        end
    end

    if CalculateCD == 1 then
        if CalculationRestrictions == nil then
            Gr = Gr + CreateSpectra(H_f, Tr_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}}) * dZ
            Gl = Gl + CreateSpectra(H_f, Tl_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}}) * dZ
        else
            Gr = Gr + CreateSpectra(H_f, Tr_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'restrictions', CalculationRestrictions}}) * dZ
            Gl = Gl + CreateSpectra(H_f, Tl_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'restrictions', CalculationRestrictions}}) * dZ
        end
    end

    if CalculateLD == 1 then
        if CalculationRestrictions == nil then
            Gein1 = Gein1 + CreateSpectra(H_f, Tein1_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}}) * dZ
            Gein2 = Gein2 + CreateSpectra(H_f, Tein2_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}}) * dZ
        else
            Gein1 = Gein1 + CreateSpectra(H_f, Tein1_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'restrictions', CalculationRestrictions}}) * dZ
            Gein2 = Gein2 + CreateSpectra(H_f, Tein2_2p_3d, Psi, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'restrictions', CalculationRestrictions}}) * dZ
        end
    end
end
io.write(string.format('===============\n'))

LorentzianFWHM_L3 = 0.4 - Gamma
LorentzianFWHM_L2 = 0.4 - Gamma
Emiddle = 862.7 - DeltaE

GaussianFWHM = 2 * math.sqrt(2 * math.log(2)) * 0.1

if CalculateIso == 1 then
    Giso = Giso / Z / 3
    Giso.Broaden(GaussianFWHM, {{Emin, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L2}, {Emax, LorentzianFWHM_L2}})
    Giso.Print({{'file', 'input' .. '_iso.spec'}})
end

if CalculateCD == 1 then
    Gr = Gr / Z
    Gl = Gl / Z
    Gcd = Gr - Gl
    Gcd.Broaden(GaussianFWHM, {{Emin, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L2}, {Emax, LorentzianFWHM_L2}})
    Gcd.Print({{'file', 'input' .. '_cd.spec'}})
    Gl.Broaden(GaussianFWHM, {{Emin, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L2}, {Emax, LorentzianFWHM_L2}})
    Gl.Print({{'file', 'input' .. '_left.spec'}})
    Gr.Broaden(GaussianFWHM, {{Emin, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L2}, {Emax, LorentzianFWHM_L2}})
    Gr.Print({{'file', 'input' .. '_right.spec'}})
end

if CalculateLD == 1 then
    Gein1 = Gein1 / Z
    Gein2 = Gein2 / Z
    Gld = Gein1 - Gein2
    Gld.Broaden(GaussianFWHM, {{Emin, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L3}, {Emiddle, LorentzianFWHM_L2}, {Emax, LorentzianFWHM_L2}})
    Gld.Print({{'file', 'input' .. '_ld.spec'}})
end

