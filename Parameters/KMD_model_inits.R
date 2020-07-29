initParms <- function(newParms = NULL) {
  parms <- c(
    bw = 0.0,
    qc = 0.0,
    vliv = 0.0,
    vpls = 0.0,
    vrbc = 0.0,
    hct = 0.42,
    vspf = 0.0,
    vrpf = 0.0,
    rvurine = 0,
    vurinec = 0.0,
    uclmet = 0.0,
    qliv = 0.0,
    qspf = 0.0,
    qrpf = 0.0,
    pliv = 0.0,
    prpf = 0.0,
    pspf = 0.0,
    vmax = 0.0,
    km = 0.0,
    vkm1 = 0.0,
    MW = 0.0,
    ivdoseum = 0,
    ivdose = 0.0,
    ivlen = 0.0,
    odoseum = 0,
    odose = 0.0,
    ka = 0.0
  )

  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
    parms[names(newParms)] <- newParms
  }

  parms <- within(as.list(parms), {
  })
  out <- .C("getParms",  as.double(parms),
            out=double(length(parms)),
            as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

Outputs <- c(
    "cpls",
    "cliv",
    "cspf",
    "crpf",
    "mbal",
    "ramet",
    "cumet",
    "vbal",
    "qbal"
)

initStates <- function(parms, newStates = NULL) {
  Y <- c(
    apls = 0.0,
    aliv = 0.0,
    aspf = 0.0,
    arpf = 0.0,
    amet = 0.0,
    aumet = 0.0,
    aiv = 0.0
  )

  if (!is.null(newStates)) {
    if (!all(names(newStates) %in% c(names(Y)))) {
      stop("illegal state variable name in newStates")
    }
    Y[names(newStates)] <- newStates
  }

.C("initState", as.double(Y));
Y
}
