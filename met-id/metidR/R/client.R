.metid_error <- function(message, status_code = NA_integer_) {
  condition <- structure(
    list(
      message = message,
      call = NULL,
      status_code = status_code
    ),
    class = c("metid_api_error", "error", "condition")
  )
  stop(condition)
}

.metid_repo_path <- function(...) {
  root <- Sys.getenv("METID_PROJECT_ROOT", unset = "")

  if (!nzchar(root)) {
    stop("Set METID_PROJECT_ROOT in .Renviron.")
  }

  file.path(root, ...)
}


.metid_request <- function(client, method, path, body = NULL,
                           retries = 3, timeout = 60) {
  url <- paste0(sub("/$", "", client$base_url), path)
  last_error <- NULL

  for (attempt in seq_len(retries)) {
    response <- tryCatch(
      httr::VERB(
        method = method,
        url = url,
        body = body,
        encode = "json",
        httr::timeout(timeout)
      ),
      error = function(e) e
    )

    if (inherits(response, "error")) {
      last_error <- response$message

      if (attempt < retries) {
        Sys.sleep(0.5)
        next
      }

      .metid_error(
        paste("Could not contact Met-ID server:", last_error)
      )
    }

    status <- httr::status_code(response)
    text <- httr::content(response, as = "text", encoding = "UTF-8")

    if (status >= 400) {
      message <- text

      payload <- tryCatch(
        jsonlite::fromJSON(text, simplifyVector = FALSE),
        error = function(e) NULL
      )

      if (is.list(payload)) {
        message <- payload$error %||% payload$message %||% text
      }

      if (status < 500 || attempt == retries) {
        .metid_error(
          paste0("Met-ID API returned HTTP ", status, ": ", message),
          status
        )
      }

      last_error <- message
      Sys.sleep(0.5)
      next
    }

    if (!nzchar(text)) {
      return(NULL)
    }

    return(jsonlite::fromJSON(text, simplifyVector = FALSE))
  }

  .metid_error(paste("Met-ID request failed:", last_error))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

.metid_table_frame <- function(table) {
  column_names <- unlist(table$col_names, use.names = FALSE)
  rows <- table$rows

  if (length(rows) == 0) {
    result <- as.data.frame(
      setNames(replicate(length(column_names), logical(0), simplify = FALSE),
               column_names),
      check.names = FALSE
    )
    return(result)
  }

  columns <- lapply(seq_along(column_names), function(j) {
    values <- lapply(rows, function(row) row[[j]])

    values <- lapply(values, function(value) {
      if (is.null(value)) NA else value
    })

    unlist(values, recursive = FALSE, use.names = FALSE)
  })

  names(columns) <- column_names

  as.data.frame(
    columns,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

metid_client <- function(base_url = "http://127.0.0.1:8787") {
  structure(
    list(
      base_url = sub("/$", "", base_url),
      process = NULL
    ),
    class = "metid_client"
  )
}

start_metid_server <- function(
    binary = NULL,
    sidecar = NULL,
    port = 8787,
    timeout = 60
) {
  # Explicit arguments take priority, then .Renviron overrides.
  server_env <- Sys.getenv("METID_SERVER_PATH", unset = "")
  if (is.null(binary) && nzchar(server_env)) {
    binary <- server_env
  }

  sidecar_env <- Sys.getenv("METID_CCS_SIDECAR_PATH", unset = "")
  if (is.null(sidecar) && nzchar(sidecar_env)) {
    sidecar <- sidecar_env
  }

  if (is.null(binary)) {
    executable <- if (.Platform$OS.type == "windows") {
      "Met-ID.exe"
    } else {
      "Met-ID"
    }

    binary <- .metid_repo_path(
      "met-id", "src-tauri", "target", "release", executable
    )
  }

  if (is.null(sidecar)) {
    pattern <- if (.Platform$OS.type == "windows") {
      "prediction_model-*.exe"
    } else {
      "prediction_model-*"
    }

    candidates <- Sys.glob(.metid_repo_path(
      "met-id", "src-tauri", "pyinstaller", "dist", pattern
    ))

    if (length(candidates) != 1) {
      stop(
        "Expected to find one CCS sidecar, but found ",
        length(candidates),
        ". Set METID_CCS_SIDECAR_PATH in .Renviron."
      )
    }

    sidecar <- candidates[[1]]
  }

  if (!file.exists(binary)) {
    stop("Server binary not found: ", binary)
  }

  if (!file.exists(sidecar)) {
    stop("CCS sidecar not found: ", sidecar)
  }

  binary <- normalizePath(binary, winslash = "/", mustWork = TRUE)
  sidecar <- normalizePath(sidecar, winslash = "/", mustWork = TRUE)

  process <- processx::process$new(
    command = binary,
    args = c("--server", "--port", as.character(port)),
    env = c(METID_CCS_SIDECAR = sidecar),
    stdout = NULL,
    stderr = NULL,
    cleanup = TRUE
  )

  client <- structure(
    list(
      base_url = paste0("http://127.0.0.1:", port),
      process = process
    ),
    class = "metid_client"
  )

  deadline <- Sys.time() + timeout

  repeat {
    if (!process$is_alive()) {
      stop("Met-ID server exited before becoming ready.")
    }

    ready <- tryCatch(
      httr::GET(
        paste0(client$base_url, "/warmup"),
        httr::timeout(1)
      ),
      error = function(e) NULL
    )

    if (!is.null(ready) && httr::status_code(ready) == 200) {
      return(client)
    }

    if (Sys.time() >= deadline) {
      process$kill()
      stop("Timed out waiting for the Met-ID server to start.")
    }

    Sys.sleep(0.1)
  }
}

stop_metid_server <- function(client) {
  if (!is.null(client$process) && client$process$is_alive()) {
    client$process$kill()
  }

  client$process <- NULL
  invisible(client)
}

metid_name <- function(client, name, retries = 3) {
  .metid_request(
    client,
    "GET",
    "/name",
    body = name,
    retries = retries
  )
}

metid_process <- function(client, items, retries = 3) {
  if (!is.numeric(items)) {
    stop("`items` must be a numeric vector.")
  }

  .metid_request(
    client,
    "POST",
    "/process",
    body = as.list(items),
    retries = retries
  )
}

metid_identify <- function(
    client,
    metabolome = "HMDB(All)",
    matrix = "FMP-10",
    typ = "All",
    adducts = "All",
    mass_error = 0,
    mz = 421.1911,
    mzwindow = 5,
    ccs = -1,
    retries = 1
) {
  if (!is.numeric(mz) || length(mz) == 0) {
    stop("`mz` must be a non-empty numeric vector.")
  }

  if (is.null(ccs)) {
    ccs <- rep(-1, length(mz))
  } else if (length(ccs) == 1 && ccs == -1) {
    ccs <- rep(-1, length(mz))
  }

  if (length(ccs) != length(mz)) {
    stop("Provide one CCS value per m/z value.")
  }

  payload <- list(
    met = metabolome,
    mat = matrix,
    typ = typ,
    adducts = adducts,
    mass_error = mass_error,
    masses = as.list(mz),
    ccss = as.list(ccs),
    mzwindow = mzwindow
  )

  result <- .metid_request(
    client,
    "POST",
    "/identify",
    body = payload,
    retries = retries
  )

  columns <- lapply(result, unlist, use.names = FALSE)
  output <- as.data.frame(
    columns,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if ("names" %in% names(output)) {
    rownames(output) <- output$names
  }

  output
}

metid_list_tables <- function(client, retries = 3) {
  result <- .metid_request(
    client,
    "GET",
    "/get_database",
    body = "tables",
    retries = retries
  )

  unlist(result$Tables, use.names = FALSE)
}

metid_get_table <- function(client, table, retries = 3) {
  result <- .metid_request(
    client,
    "GET",
    "/get_database",
    body = table,
    retries = retries
  )

  if (is.null(result$Table)) {
    stop("The server did not return a database table.")
  }

  .metid_table_frame(result$Table)
}

metid_list_metabolites <- function(client, retries = 3) {
  metid_get_table(client, "metabolites", retries)
}

metid_list_matrices <- function(client, retries = 3) {
  metid_get_table(client, "matrices", retries)
}

metid_list_functional_groups <- function(client, retries = 3) {
  metid_get_table(client, "functional_groups", retries)
}

metid_list_adducts <- function(client, retries = 3) {
  metid_get_table(client, "adducts", retries)
}

metid_search_by_name <- function(client, name, ccs = FALSE, retries = 3) {
  result <- .metid_request(
    client,
    "GET",
    "/search_by_name",
    body = list(name = name, ccs = ccs),
    retries = retries
  )

  if (is.null(result$Table)) {
    stop("The server did not return search results.")
  }

  output <- .metid_table_frame(result$Table)

  if ("name" %in% names(output)) {
    rownames(output) <- output$name
  }

  output
}

metid_read_imzml <- function(client, path, n, tolerance, retries = 3) {
  if (!file.exists(path)) {
    stop("imzML file not found: ", path)
  }

  result <- .metid_request(
    client,
    "POST",
    "/read_imzml",
    body = list(
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      n = n,
      tolerance = tolerance
    ),
    retries = retries
  )

  data.frame(
    mz = unlist(result$mz, use.names = FALSE),
    intensity = unlist(result$intensity, use.names = FALSE)
  )
}

metid_predict_ccs <- function(
    client,
    smiles,
    adduct = NULL,
    retries = 3
) {
  if (!is.character(smiles) || length(smiles) == 0 ||
      anyNA(smiles) || any(!nzchar(trimws(smiles)))) {
    stop("Provide one or more non-empty SMILES strings.")
  }

  payload <- list(smiles = as.list(smiles))

  if (!is.null(adduct)) {
    payload$adducts <- as.list(adduct)
  }

  result <- .metid_request(
    client,
    "POST",
    "/predict_ccs",
    body = payload,
    retries = retries
  )

  if (length(result) == 0) {
    return(data.frame())
  }

  rows <- lapply(result, function(row) {
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  })

  do.call(rbind, rows)
}