check_prev_out <- function(prev_out, check_names) {
  for (check_name in check_names) {
    if (check_name %in% names(prev_out)) next()

    stop(check_name, " is missing from prev_out but is required.")
  }
}

check_input <- function(input) {

  # check that metadata items length is same as number of samples
  metadata <- input$metadata

  if (length(metadata)) {
    nsamples <- length(input$sampleNames)
    nmeta <- sapply(metadata, length)
    if (!all(nmeta == nsamples)) stop("Sample number differs from metadata length.")
  }
}

get_color_pool <- function() {
  c(
    "#77aadd", "#ee8866", "#eedd88", "#ffaabb", "#99ddff", "#44bb99",
    "#bbcc33", "#999900", "#bc9cc9", "#295a8b", "#b14d2c", "#c2aa36",
    "#d4788a", "#3e8ab1", "#158364", "#737a38", "#4e4e00", "#835794",
    "#edf761", "#4d8cc4", "#c7221f", "#6eb288", "#74008b", "#b3fd33",
    "#d7ae3e", "#526fbd", "#5da8a2", "#9bfdf8", "#a473b0", "#e59637",
    "#c2c2c2", "#521913", "#e04e29", "#529cb4", "#8b201a", "#8fbc65",
    "#fdb462", "#bebada", "#ff5a48", "#fccde5", "#bc80bd", "#ccebc5",
    "#6a3d9a", "#9b7db7", "#cafffe", "#8dd3c7", "#1f78b4", "#0178ff",
    "#ff6db6", "#fff5c3", "#b8787a", "#33a02c", "#52bff8", "#e89976",
    "#ff7f00", "#ffecd1", "#19c2c1", "#ffff50", "#bda985", "#5bf1fc",
    "#98b299", "#d1f8b0", "#51bcae", "#ebffc5", "#c1ffd2", "#e0ffb3",
    "#79baff", "#fdbf6f", "#ff9587", "#a6cee3", "#ffde5d", "#a8d5ff",
    "#ffb05d", "#fff48c", "#9183a7", "#469966", "#e31a1c", "#c4c7ff",
    "#c3781d", "#7eeece", "#d9d9d9", "#f5aae2", "#c0fdbd", "#5ab9c2",
    "#feff84", "#bcab6a", "#60afff", "#dddddd", "#809055", "#b15928",
    "#7b7fe8", "#cdf1e1", "#ffdafd", "#5b9667", "#ecec75", "#6f9528",
    "#e1ff82", "#03dec8", "#d2c233", "#e38423", "#94ffc9", "#01c2e2",
    "#b7ffcf", "#639650", "#ffbe59", "#fff0c6", "#abb094", "#c7752e",
    "#ffdd6f", "#debc76", "#03d4a0", "#6bb4da", "#01bba7", "#ffa86c",
    "#cd63b3", "#7e9061", "#ff9b54", "#ff6b76", "#90e065", "#ffaaa6",
    "#e2ffd0", "#4d9b32", "#ffffaf", "#7d86c1", "#01be8b", "#1a8bf1",
    "#8d83b3", "#ffcab6", "#af6ddb", "#84e672", "#329a7a", "#02fdfc",
    "#6ebb97", "#91b585", "#b77e40", "#2dbfff", "#77ffbd", "#ffb5fd",
    "#cab2d6", "#d76a39", "#8579f3", "#a3bdff", "#f5c23e", "#fff1b5",
    "#f6ffd0", "#aadaff", "#ac831b", "#c4764b", "#ffdeae", "#a8b624",
    "#c87505", "#f16f3c", "#ffe07e", "#b67f07", "#8b8e3b", "#5b99ff",
    "#da674b", "#9fc83f", "#bec838", "#58957d", "#688f9e", "#e35cc2",
    "#01edc1", "#7289bc", "#ff9eb0", "#ffdec2", "#ffd2d3", "#f8c8ff",
    "#02a287", "#be7672", "#02fcba", "#6089d3", "#b58500", "#01d68a",
    "#a6b9ff", "#4ef6ff", "#b0d64c", "#65b237", "#d4ffdd", "#cb9fff",
    "#ffbfe6", "#7aff9b", "#789424", "#84ffb2", "#66973c", "#aa833f",
    "#ffdff7", "#ff6bb5", "#d2c2a8", "#ca7423", "#ff9064", "#b0d2ff",
    "#caddc3", "#ff9174", "#00a96f", "#d3b1a4", "#d4d243", "#1ba33e",
    "#45ffe5", "#dc6178", "#23989e", "#edff81", "#d166d6", "#ffa8b2",
    "#ff999d", "#00dfa5", "#4f9a4c", "#fff26d", "#f6ffbb", "#00ada6",
    "#4dda79", "#a478c0", "#4589e9", "#ebb6ff", "#7ca4ff", "#a2ffb1",
    "#819900", "#bfffe4", "#918d14", "#ffd99f", "#8394ff", "#01c7d5",
    "#f1da4d", "#01eda8", "#d16d4f", "#01f3bd", "#b97a66", "#b2df8a",
    "#bd7958", "#a98074", "#a1ffdc", "#72cfff", "#0197cb", "#01b0f3",
    "#6e8dab", "#36e5ff", "#988868", "#3d8fcd", "#fd508c", "#cdff8d",
    "#a490ff", "#53e1ff", "#fa9dff", "#359d43", "#00afb4", "#bca80f",
    "#02c5c9", "#eeffa7", "#f6ff98", "#ffaa4c", "#d98816", "#ff84ec",
    "#abc335", "#dfe1ff", "#b07995", "#02c8f9", "#3591c5", "#609e18",
    "#029fc8", "#70de6f", "#01cfa8", "#c168c0", "#7883d6", "#cdffa4",
    "#e4539a", "#01a1c0", "#02b68d", "#ffa779", "#059e52", "#ded5ff",
    "#85ffa4", "#aba605", "#679282", "#019dae", "#ce703e", "#01c271",
    "#d7ffc0", "#ffb675", "#46968d", "#6bc752", "#fffe71", "#ff736b",
    "#649196", "#848d6c", "#a17bf5", "#99fff9", "#5d86e9", "#977bce",
    "#9e883d", "#dc96ff", "#aef473", "#88fff1", "#d4ff99", "#9cff9d",
    "#01e99a", "#da7df1", "#a47cab", "#e99428", "#93d453", "#53ffb5",
    "#98f9ff", "#f0b334", "#378ddc", "#ffad96", "#d26b61", "#c8baff",
    "#819136", "#7befff", "#ffb3e8", "#d0ae1f", "#d6f467", "#ffe395",
    "#ccfff3", "#ffb98d", "#b5b2ff", "#ceffd3", "#8ec13b", "#c4fd77",
    "#02b6d5", "#00fef2", "#d6b1ff", "#d06c6b", "#ceffb1", "#00aa4f",
    "#9a9400", "#01f7ce", "#98e5ff", "#fa8a3b", "#97c4ff", "#ef564b",
    "#6a9375", "#4093af", "#af7f5c", "#3b9a72", "#ba72a1", "#ff625c",
    "#c49d04", "#9e8920", "#ff6071", "#ff8cc6", "#02a574", "#bf7488",
    "#ff9a41", "#7aafff", "#f6f0d2", "#c2fcff", "#33d5ff", "#af72c2",
    "#4c9958", "#4eba4d", "#8dd7ff", "#ed7434", "#ffcba4", "#7dffd4",
    "#ffeb5d", "#ffeea6", "#c07e00", "#ffc682", "#9f8654", "#ff8db5",
    "#7ea109", "#46e98e", "#00d8df", "#1fecff", "#918c42", "#ffcf74",
    "#3cc360", "#a0f87e", "#ff638d", "#b0efff", "#bf7940", "#a08094",
    "#fb9a99", "#ffabcb", "#009f91", "#d1678d", "#a8838a", "#ff9ee7",
    "#f977e0", "#f46149", "#00edd6", "#7e913f", "#ffc0cb"
  )
}
