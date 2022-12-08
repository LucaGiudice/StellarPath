#' Initialize StellarPath required python packages
#'
#' Initialize StellarPath required python packages to call from R.
#' Credits for this function go to spacyr R package: https://github.com/quanteda/spacyr/tree/master/R
#' Respect for these fellow R developers who developed a way to automatically create and detect a conda enviroment linked to
#' an R package. Hope that Reticulate will integrate this function in their standard package for the future.
#'
#' @return NULL
#' @param python_executable the full path to the Python executable, for which
#'   StellarPath required python packages is installed.
#' @param ask logical; if \code{FALSE}, use the first StellarPath required python packages installation found;
#'   if \code{TRUE}, list available StellarPath required python packages installations and prompt the user for
#'   which to use. If another (e.g. \code{python_executable}) is set, then this
#'   value will always be treated as \code{FALSE}.
#' @param virtualenv set a path to the Python virtual environment with StellarPath required python packages
#'   installed Example: \code{virtualenv = "~/myenv"}
#' @param condaenv set a path to the anaconda virtual environment with StellarPath required python packages
#'   installed Example: \code{condalenv = "myenv"}
#' @param check_env logical; check whether conda/virtual environment generated
#'   by \code{SP_install()} exists
#' @param refresh_settings logical; if \code{TRUE}, StellarPath will ignore the saved
#'   settings in the profile and initiate a search of new settings.
#' @param save_profile logical; if \code{TRUE}, the current StellarPath required python packages setting will
#'   be saved for the future use.
#' @param prompt logical; asking whether user wants to set the environment as default.
#' @export
SP_initialize <- function(python_executable = NULL,
                               virtualenv = NULL,
                               condaenv = "SP_condaenv",
                               ask = FALSE,
                               refresh_settings = FALSE,
                               save_profile = FALSE,
                               check_env = TRUE,
                               prompt = TRUE) {
  set_SP_python_option(
    python_executable,
    virtualenv,
    condaenv,
    check_env,
    refresh_settings,
    ask
  )

  ## check settings and start reticulate python
  settings <- check_SP_python_options()
  if (!is.null(settings)) {
    if (settings$key == "SP_python_executable") {
      reticulate::use_python(settings$val, required = TRUE)
    } else if (settings$key == "SP_virtualenv") {
      reticulate::use_virtualenv(settings$val, required = TRUE)
    } else if (settings$key == "SP_condaenv") {
      reticulate::use_condaenv(settings$val, required = TRUE)
    }
  }

  # Importing this here may start importing necessary packages
  reticulate::source_python(system.file("python",
    "sage_classification_individual.py",
    package = "StellarPath",
    mustWork = TRUE
  ))

  message(colourise(
    "\nSuccessfully initialized StellarPath required python packages.\n",
    fg = "green", bg = NULL
  ))
  settings <- check_SP_python_options()

  settings_SP <- paste('Python options: \n type = "', settings$key,
    '", \n name = "', settings$val, '".',
    sep = ""
  )

  message(colourise(settings_SP,
    fg = "blue", bg = NULL
  ))


  options("SP_initialized" = TRUE)

  if (save_profile == TRUE) {
    save_SP_options(settings$key, settings$val, prompt = prompt)
  }

}

#' Find StellarPath required python packages
#'
#' Locate the user's version of Python for which StellarPath required python packages are installed.
#' @return SP_python
#' @export
#' @param ask logical; if \code{FALSE}, use the first StellarPath required python packages installation found;
#'   if \code{TRUE}, list available StellarPath required python packages installations and prompt the user
#'   for which to use. If another (e.g. \code{python_executable}) is set, then
#'   this value will always be treated as \code{FALSE}.
#'
#' @keywords internal
find_SP <- function(ask) {
  SP_found <- `:=` <- NA
  SP_python <- NULL
  options(warn = -1)
  py_execs <- if (is_windows()) {
    system2("where", "python", stdout = TRUE)
  } else if (is_osx() && file.exists("~/.bash_profile")) {
    c(
      system2("source", "~/.bash_profile; which -a python", stdout = TRUE),
      system2("source", "~/.bash_profile; which -a python3", stdout = TRUE)
    )
  } else {
    c(
      system2("which", "-a python", stdout = TRUE),
      system2("which", "-a python3", stdout = TRUE)
    )
  }
  py_execs <- unique(py_execs)
  options(warn = 0)

  if (length(py_execs) == 0 | grepl("not find", py_execs[1])[1]) {
    return(NA)
  }

  df_python_check <- tibble::tibble(py_execs, SP_found = 0)
  for (i in seq_len(nrow(df_python_check))) {
    py_exec <- df_python_check[i, ]
    sys_message <- check_SP_model(py_exec) #
    if (sys_message == "OK") {
      df_python_check[i, SP_found := 1]
    }
  }

  if (df_python_check[, sum(SP_found)] == 0) {
    return(NULL)
  } else if (df_python_check[, sum(SP_found)] == 1) {
    SP_python <- df_python_check[SP_found == 1, py_execs]
    message("StellarPath: ", ") is installed in ", SP_python)
  } else if (ask == FALSE) {
    SP_python <- df_python_check[SP_found == 1, py_execs][1]
    message("StellarPath: is installed in more than one python")
    message("StellarPath will use ", SP_python, " (because ask = FALSE)")
  } else {
    SP_pythons <- df_python_check[SP_found == 1, py_execs]
    message("StellarPath is installed in more than one python")
    number <- utils::menu(SP_pythons, title = "Please select python:")
    if (number == 0) {
      stop("Initialization was canceled by user", call. = FALSE)
    }
    SP_python <- SP_pythons[number]
    message("StellarPath will use: ", SP_python)
  }
  return(SP_python)
}


#' Find StellarPath required python pacakges env
#'
#' check whether conda/virtual environment for StellarPath required python pacakges exists
#' @export
#'
#' @keywords internal
find_SP_env <- function() {
  if (is.null(tryCatch(reticulate::conda_binary("auto"), error = function(e) NULL))) {
    return(FALSE)
  }
  found <- if ("SP_condaenv" %in% reticulate::conda_list(conda = "auto")$name) {
    TRUE
  } else if (file.exists(file.path("~/.virtualenvs", "SP_virtualenv", "bin", "activate"))) {
    TRUE
  } else {
    FALSE
  }
  return(found)
}


check_SP_model <- function(py_exec) { ### , model
  options(warn = -1)
  py_exist <- if (is_windows()) {
    if (py_exec %in% system2("where", "python", stdout = TRUE)) {
      py_exec
    } else {
      NULL
    }
  } else {
    system2("which", py_exec, stdout = TRUE)
  }

  if (length(py_exist) == 0) {
    stop(py_exec, " is not a python executable")
  }
  tryCatch({
    sys_message <- "see error in SP_initialize row 235"
    # system2(py_exec, c(sprintf("-c \"import texrpp; StellarPath.load('%s'); print('OK')\"", model)),
    #        stderr = TRUE, stdout = TRUE)
  })
  options(warn = 0)
  return(paste(sys_message, collapse = " "))
}


set_SP_python_option <- function(python_executable = NULL,
                                      virtualenv = NULL,
                                      condaenv = NULL,
                                      check_env = TRUE,
                                      refresh_settings = FALSE,
                                      ask = NULL) {
  if (refresh_settings) clear_SP_options()

  if (!is.null(check_SP_python_options())) {
    settings <- check_SP_python_options()

    message_SP1 <- paste("StellarPath python option is already set, StellarPath will use: ",
      sub("SP_", "", settings$key), ' = "', settings$val, '"',
      sep = ""
    )

    message(colourise(message_SP1,
      fg = "blue", bg = NULL
    ))
  }
  # a user can specify only one
  else if (sum(!is.null(c(python_executable, virtualenv, condaenv))) > 1) {
    stop(paste(
      "Too many python environments are specified, please select only one",
      "from python_executable, virtualenv, and condaenv"
    ))
  }
  # give warning when nothing is specified
  else if (sum(!is.null(c(python_executable, virtualenv, condaenv))) == 1) {
    if (!is.null(python_executable)) {
      if (check_SP_model(python_executable) != "OK") {
        stop("StellarPath required python packages ", " are not installed in ", python_executable)
      }
      clear_SP_options()
      options(SP_python_executable = python_executable)
    } else if (!is.null(virtualenv)) {
      clear_SP_options()
      options(SP_virtualenv = virtualenv)
    } else if (!is.null(condaenv)) {
      clear_SP_options()
      options(SP_condaenv = condaenv)
    }
  } else if (check_env &&
    !(is.null(tryCatch(reticulate::conda_binary("auto"), error = function(e) NULL))) &&
    "SP_condaenv" %in% reticulate::conda_list(conda = "auto")$name) {
    message(colourise(
      "Found 'SP_condaenv'. StellarPath will use this environment \n",
      fg = "green", bg = NULL
    ))
    clear_SP_options()
    options(SP_condaenv = "SP_condaenv")
  } else if (check_env && file.exists(file.path("~/.virtualenvs", virtualenv, "bin", "activate"))) {
    message(colourise(
      "Found your specified virtual environment. StellarPath will use this environment \n",
      fg = "green", bg = NULL
    )) # OK: original: Found 'SP_virtualenv'. StellarPath will use this environment"
    clear_SP_options()
    options(SP_virtualenv = file.path("~/.virtualenvs/", virtualenv))
  } else {
    message("Finding a python executable with StellarPath required python pakages installed...")
    SP_python <- find_SP(ask = ask) # model,
    if (is.null(SP_python)) {
      stop("StellarPath required python packages ", " are not installed in any of python executables.") #  model,
    } else if (is.na(SP_python)) {
      stop("No python was found on system PATH")
    } else {
      options(SP_python_executable = SP_python)
    }
  }
  return(NULL)
}


clear_SP_options <- function() {
  options(SP_python_executable = NULL)
  options(SP_condaenv = NULL)
  options(SP_virtualenv = NULL)
}

check_SP_python_options <- function() {
  settings <- NULL
  for (k in c(
    "SP_python_executable",
    "SP_condaenv",
    "SP_virtualenv"
  )) {
    if (!is.null(getOption(k))) {
      settings$key <- k
      settings$val <- getOption(k)
    }
  }
  return(settings)
}

save_SP_options <- function(key, val, prompt = TRUE) {
  prof_file <- "~/.Rprofile"
  if (!is.null(getOption("SP_prompt"))) prompt <- getOption("SP_prompt")

  ans <- if (prompt) {
    utils::menu(c("No", "Yes"),
      title = sprintf('Do you want to set the option, \'%s = "%s"\' , as a default (y|[n])? ', key, val)
    )
  } else {
    2
  }
  if (ans == 2) {
    rprofile <- if (file.exists(prof_file)) readLines(prof_file) else NULL
    rprofile <- grep("options\\(\\s*SP_.+\\)", rprofile, value = TRUE, invert = TRUE)
    rprofile <- c(rprofile, sprintf('options(%s = "%s")', key, val))
    write(rprofile, file = prof_file)
    message(colourise(
      "The option was saved. The option will be used in SP_initialize() in future \n",
      fg = "green", bg = NULL
    ))
  } else {
    message("The option was not saved (user cancelled)")
  }
}

#' Test StellarPath python script
#'
#' Test StellarPath python script
#' @return NULL
#' @importFrom reticulate source_python
#' @export
test_SP_py = function(){
  # Importing this here may start importing necessary packages
  reticulate::source_python(system.file("python",
                                        "sage_classification_individual.py",
                                        package = "StellarPath",
                                        mustWork = TRUE
  ))
  message(colourise(
    "Reticulate enviroment and the StellarPath python module work perfectly",
    fg = "green", bg = NULL
  ))
}
