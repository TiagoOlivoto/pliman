# Progress bar
# used in metan R package
# https://github.com/TiagoOlivoto/metan/blob/master/R/utils_progress.R
sec_to_hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":"
  )
}
progress <- function(min = 0,
                     max = 100,
                     leftd = "|",
                     rightd = "|",
                     char = "=",
                     style = 2,
                     width = getOption("width"),
                     time = Sys.time()){
  # Adapted from https://stackoverflow.com/a/26920123/15245107
  return(list(min = min,
              max = max,
              leftd = leftd,
              rightd = rightd,
              char = char,
              style = style,
              width = width,
              time = time))
}
run_progress <- function(pb,
                         actual,
                         text = "",
                         digits = 0,
                         sleep = 0){
  Sys.sleep(sleep)
  elapsed <- sec_to_hms(as.numeric(difftime(Sys.time(), pb$time, units = "secs")))
  temp <- switch(
    pb$style,
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd),
         text = paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""))),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 6,
         text =  paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""), '% s%%')),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 9,
         text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), elapsed)),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 15,
         text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), '% s%%', elapsed))
  )
  step <- round(actual / pb$max * (pb$width - temp$extra))
    cat(sprintf(temp$text,
                strrep(pb$char, step),
                strrep(' ', pb$width - step - temp$extra),
                round(actual / pb$max * 100, digits = digits)), "\r")
  if(actual == pb$max){
    cat("\n")
  }
}

check_names_dir <- function(name, names_dir, dir){
  if(!name %in% names_dir){
    stop(paste("'", name, "' not found in '",
               paste(getwd(), sub(".", "", dir), sep = ""), "'", sep = ""),
         call. = FALSE)
  }
}

image_correct <- function(image, perc){
  t <- image
  n=round(perc*min(c(ncol(t),nrow(t))),0)
  p1=function(t){
    t2=t
    for( i in 2:(nrow(t)-n-1)){
      for(j in 1:ncol(t)){
        if(t[i,j]==1){
          if(t[i-1,j]==0){
            a=0
            while(a<n){
              a=a+1

              if(sum(t[i:(i+a),j]==1)<a){t2[i:(i+a),j]=0;a=n}
            }
          }
        }
      }
    }
    return(t2)
  }
  Pp=p1(t)
  Pp=p1(t(Pp))
  return(t(Pp))
}

