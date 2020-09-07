import numpy as _np


# TODO: Documentation!
def Get_CLVD_DC(Full_Moment):
    M = Full_Moment
    # isotropic part
    M_iso = _np.diag(
        _np.array([1.0 / 3 * _np.trace(M), 1.0 / 3 * _np.trace(M), 1.0 / 3 * _np.trace(M)])
    )

    M0_iso = abs(1.0 / 3 * _np.trace(M))

    # deviatoric part
    M_devi = M - M_iso

    isotropic = M_iso
    deviatoric = M_devi

    # eigenvalues and -vectors
    eigenwtot, eigenvtot = _np.linalg.eig(M_devi)

    # eigenvalues and -vectors of the deviatoric part
    eigenw1, eigenv1 = _np.linalg.eig(M_devi)

    # eigenvalues in ascending order:
    eigenw = _np.real(_np.take(eigenw1, _np.argsort(abs(eigenwtot))))
    eigenv = _np.real(_np.take(eigenv1, _np.argsort(abs(eigenwtot)), 1))

    # eigenvalues in ascending order in absolute value!!:
    eigenw_devi = _np.real(_np.take(eigenw1, _np.argsort(abs(eigenw1))))
    # eigenv_devi = _np.real(_np.take(eigenv1, _np.argsort(abs(eigenw1)), 1))

    M0_devi = max(abs(eigenw_devi))

    # named according to Jost & Herrmann:
    # a1 = eigenv[:, 0]
    a2 = eigenv[:, 1]
    a3 = eigenv[:, 2]

    # if only isotropic part exists:
    epsilon = 1e-13
    if M0_devi < epsilon:
        F = 0.5
    else:
        F = -eigenw_devi[0] / eigenw_devi[2]

    M_DC = _np.matrix(_np.zeros((9), float)).reshape(3, 3)
    M_CLVD = _np.matrix(_np.zeros((9), float)).reshape(3, 3)

    M_DC = eigenw[2] * (1 - 2 * F) * (_np.outer(a3, a3) - _np.outer(a2, a2))
    M_CLVD = M_devi - M_DC

    # from obspy.imaging.beachball import beachball
    return M_CLVD, M_DC, F


def TDL(AN, BN):
    XN = AN[0]
    YN = AN[1]
    ZN = AN[2]
    XE = BN[0]
    YE = BN[1]
    ZE = BN[2]
    AAA = 1.0e-06
    CON = 57.2957795
    if abs(ZN) < AAA:
        FD = 90.0
        AXN = abs(XN)
        if AXN > 1.0:
            AXN = 1.0
        FT = _np.arcsin(AXN) * CON
        ST = -XN
        CT = YN
        if ST >= 0.0 and CT < 0:
            FT = 180.0 - FT
        if ST < 0.0 and CT <= 0:
            FT = 180.0 + FT
        if ST < 0.0 and CT > 0:
            FT = 360.0 - FT
        FL = _np.arcsin(abs(ZE)) * CON
        SL = -ZE
        if abs(XN) < AAA:
            CL = XE / YN
    else:
        if -ZN > 1.0:
            ZN = -1.0
        FDH = _np.arccos(-ZN)
        FD = FDH * CON
        SD = _np.sin(FDH)
        if SD == 0:
            raise ValueError("Return function...")
            # return FT,FD,FL
        ST = -XN / SD
        CT = YN / SD
        SX = abs(ST)
        if SX > 1.0:
            SX = 1.0
        FT = _np.arcsin(SX) * CON
        if ST >= 0.0 and CT < 0:
            FT = 180.0 - FT
        if ST < 0.0 and CT <= 0:
            FT = 180.0 + FT
        if ST < 0.0 and CT > 0:
            FT = 360.0 - FT
        SL = -ZE / SD
        SX = abs(SL)
        if SX > 1.0:
            SX = 1.0
        FL = _np.arcsin(SX) * CON
        if ST == 0:
            CL = XE / CT
        else:
            XXX = YN * ZN * ZE / SD / SD + YE
            CL = -SD * XXX / XN
            if CT == 0:
                CL = YE / ST

        if SL >= 0.0 and CL < 0:
            FL = 180.0 - FL
        if SL < 0.0 and CL <= 0:
            FL = FL - 180.0
        if SL < 0.0 and CL > 0:
            FL = -FL
    return FT, FD, FL


def GET_sdr_from_mij(mxx, myy, mzz, mxy, mxz, myz):
    # M = _np.array([[mxx, mxy, mxz], [mxy, myy, myz], [mxz, myz, mzz]])
    M = _np.array([[mzz, mxz, myz], [mxz, mxx, mxy], [myz, mxy, myy]])

    eigenValues, eigenVectors = _np.linalg.eig(M)

    idx = _np.flip(eigenValues.argsort()[::-1])
    D = eigenValues[idx]
    V = eigenVectors[:, idx]

    D_new = _np.array([D[2], D[0], D[1]])
    V[1:2, 0:2] = -V[1:2, 0:2]
    V_new = _np.array(
        [[V[1, 2], V[1, 0], V[1, 1]], [V[2, 2], V[2, 0], V[2, 1]], [V[0, 2], V[0, 0], V[0, 1]]]
    )

    Imin = _np.argmin(D_new)
    Imax = _np.argmax(D_new)

    AE = (V_new[:, Imax] + V_new[:, Imin]) / _np.sqrt(2.0)
    AN = (V_new[:, Imax] - V_new[:, Imin]) / _np.sqrt(2.0)
    AER = _np.sqrt(AE[0] ** 2 + AE[1] ** 2 + AE[2] ** 2)
    ANR = _np.sqrt(AN[0] ** 2 + AN[1] ** 2 + AN[2] ** 2)
    AE = AE / AER
    AN = AN / ANR

    if AN[2] <= 0.0:
        AN1 = AN
        AE1 = AE
    else:
        AN1 = -AN
        AE1 = -AE
    ft, fd, fl = TDL(AN1, AE1)
    strike = 360 - ft
    dip = fd
    rake = 180 - fl
    return strike, dip, rake

