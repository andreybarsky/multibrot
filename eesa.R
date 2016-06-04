# Exploratory Event Sequence Analysis, by Andrey Barsky
# Version 0.8 - full support for vector projected critical value determination

# By default, input for this program should be a text, csv or excel file in the format:
# a b c a b c
# a b b c a b b
# etc.
# But different ways of separating events from each other, and sequences from each other,
# can be set with the eventsep and seqsep arguments.
# Excel events can be separated by cells or by spaces, and the program should detect that.

# Any way of specifying events will do - letters, words, numbers or symbols, and however many,
# but should not be the same as either separator argument (for obvious reasons). 
# The @ sign and the dash - are default operators in the recombine function.

# You can also pass different kinds of data as input using the datatype argument, below.

# datatype is an argument specifying the structure of the data you're passing to the eesa function.
# the default is 'filetype' or 'f', treating the input as a file on your machine and reads it as a table
# (it should handle text, csv and excel files in sensible formats)
# if it is 'vector' or 'v', it treats it as a character vector of event sequences
# if it is 'string' or 's', it treats it as an atomic character vector separated by linebreaks

# terminals is an optional argument which appends an implicit START and END event to every row of events.
# It is recommended you set it to TRUE, unless your data already has explicit start and end events, 
# or you don't wish to have them for some other reason.

# critval specifies the critical value of standardised normal residual that marks particular 
# transitions as notable. There are a few ways to do this: 'projected', 'derivative' or 'chisq'.
# By default, 'projected' is calculated as the SNR with the maximum projected vector to the line connected
# the highest and lowest SNRs (above 1). This tends to find the 'elbow' of the SNR scree plot but isn't perfect,
# and is sensitive to the length of the plot's tail.
# 'derivative' calculates an estimate of the second derivative of the smoothed (LOESS) SNR curve using central tendency, and
# looks for the global maximum. This is insensitive to the plot's tail, but kind of abstract as it smooths out the distribution.
# 'chisq' calculates the critical value as sqrt(chisquare / number of cells), a cutoff recommended by Colgan & Smith (1978, p.160).
# you can also set the parameter to any valid float, though it's preferred to have a non-arbitrary cutoff.


# as 

# min.obs specifies the minimum number of times a transition must be observed
# for it to be considered notable. Ideally you want your notable transitions to be
# in the double digits but that is not possible for all data.

# collapse is an optional argument that allows for collapsing low-frequency events into a single dummy event.
# This can be helpful, since low frequency events are not always interesting but the parts
# of the sequence in which they occur might have some significance.
# By default it is set to NULL, which does no collapsing, but you can set it to 'default'
# to collapse events that occur at least 1 stdev less frequently than the mean, or
# set it to any integer to treat that as the minimum observed frequency cutoff.



# In the transition matrices, the Y axis lists antecedents. X axis lists sequiturs.

# Output is an 'eesa' object with the following elements:
# $input: How we parsed the data provided (before collapsing, etc.)
# $events: The unique event codes detected, and their raw frequencies in the input data
# $lfe: The unique event codes identified as 'low frequency', or NULL if no collapsing was done
# $observed: Matrix of observed transitional frequencies
# $expected: Matrix of expected transitional frequencies
# $residuals: Standardised normal residuals of the observed matrix
# $chisq: Chi squared statistic for the observed matrix
# $critical: Calculated critical value for standardised normal residuals to qualify as 'notable'
# $hft: Summary frame of notably high-frequency transitions (this is mainly what you want)
# $lft: Summary frame of notably low-frequency transitions

# You can then run the eesa.recombine() function on the eesa output object to
# generate a new input of higher-order events created by combining those that
# occur most frequently, and running the default eesa script again on the output
# from that function.


# enter the filepath of your own data here: 

arw_data = '/Users/lpxmac/Dropbox/Python/arw_data.txt'


eesa <- function(x, ...) UseMethod('eesa')

eesa.default <- function(x, datatype='filename', 
                         terminals=T, 
                         eventsep=' ', 
                         seqsep='\n', 
                         critval='projected', 
                         min.obs=2, 
                         collapse='default')
{

  ##### different ways of parsing the input:
  
  if(datatype == 'filename' | datatype == 'file' | datatype == 'f')
  {
    if(any(grep(".csv",x))) # if the filename contains .csv it's probably a csv
    {
      seqframe = read.csv(x, stringsAsFactors=F, header=F)
      seqvec = gsub(" NA", "", apply(seqframe, 1, paste, collapse=eventsep))
    }
    else if(any(grep(".xls",x)))  # if it contains .xls it's probably an excel file (but is it split by cell, or all in one col?)
    {
      if(!any(grep('xlsx', installed.packages()[,1]))) # if the xlsx package isn't installed, we install it
      {
        cat('Performing first time setup for reading excel files (you should only see this once)')
        install.packages("xlsx")
      }
      
      require("xlsx")
      seqframe = read.xlsx(x, 1, header=F)
      if(ncol(seqframe)>1)
      {
        seqvec = gsub(" NA", "", apply(seqframe, 1, paste, collapse=eventsep)) # collapse the df and get rid of NAs
      }
      else
      {
        seqvec = seqframe[,1]
      }
    }
    else # assume raw text, space separated
    {
      seqvec = read.table(x, header=F, sep=seqsep, stringsAsFactors=F)[,1]  
    }
    
    
  }
  else if(datatype == 'string' || datatype == 'str' | datatype == 's')
  {
    seqvec = strsplit(x, seqsep)[[1]]
  }
  else if(datatype == 'vector' || datatype == 'v')
  {
    seqvec = x
  }


  else
  {
    warning("Invalid datatype argument")
  }
  
  ##### this stuff only gets used if collapsing the low frequency events:
  
  countevents <- function(seqvec) # outputs table of raw event frequencies, for collapsing
  {
    rawstring = paste(seqvec, collapse=eventsep)
    return(table(strsplit(rawstring, eventsep)))
  }
  
  freq_cutoff <- function(eventcounts) # takes table of eventcounts, returns default cutoff value (mean - sd)
  {
    return((mean(eventcounts) - (sd(eventcounts))))
  }
  
  lfe <- function(seqvec, cutoff) # outputs vector of events that are considered low frequency
  {
    eventcounts = countevents(seqvec)
    lowfreqevents = c()
    
    for(e in 1:length(eventcounts))
    {
      if(eventcounts[e] < cutoff)
      {
        lowfreqevents = append(lowfreqevents, names(eventcounts)[e])
      }
    }
    return(lowfreqevents)
  }
  
  if(!is.null(collapse))
  {
    if(collapse == 'default') # determine cutoff value for collapsing
    {
      cutoff = freq_cutoff(countevents(seqvec))
    }
    else if(is.numeric(collapse))
    {
      cutoff = collapse
    }
    lowfreqs = lfe(seqvec, cutoff)
  }
  else
  {
    cutoff=NULL
  }
  
  ##### core functions:
  
  gettypes <- function(seqvec, cutoff_arg = cutoff, eventsep_arg = eventsep) # converts vector of sequence data into vector of discrete event types
  {
    alphabet = c()
    for(i in 1:length(seqvec))
    {
      linetypes = strsplit(seqvec[i], eventsep_arg) # break up each line
      alphabet = append(alphabet, linetypes[[1]])
    }
    eventtypes = levels(factor(alphabet)) # get the unique events
    
    if(!is.null(cutoff)) # if we are collapsing, remove the low freq events, and add the LFE event
    {
      eventtypes[!eventtypes %in% lowfreqs] # remove low freq events
      eventtypes = append(eventtypes, 'LFE') # add LFE event
    }
    
    if(length(eventtypes) == length(seqvec))
    {
      warning("Event parsing failed. Event separator may be incorrect")
      cat("Here is what was parsed:\n")
      cat(seqvec)
    }
    
    if(terminals) # if we are adding terminal events, add these to the list
    {
      return(append(eventtypes, c('START', 'END')))
    }
    else
    {
      return(eventtypes)
    }
  }
  
  
  getpairs <- function(seqvec, cutoff_arg = cutoff, eventsep_arg=eventsep) # takes sequence data vector, outputs data frame of event pairs
  {
    eventpairs = data.frame(matrix(ncol=2, nrow=1)) # a growing df with 2 cols
    i = 1
    
    for(line in 1:length(seqvec))
    {
      lineevents = strsplit(seqvec[line], eventsep_arg)[[1]]
      lineevents = lineevents[lineevents!=""] # get rid of whitespace
      
      if(length(lineevents) > 1) # this conditional fixes the case where a line is a single event
      {
        #for(e in lineevents) # fixes errors with doubled separators
        #{
        #  if(e=="")
        #  {lineevents = lineevents[!lineevents %in% e]}
        #}
        
        if(terminals) # if we are using terminals, append START here
        {
          eventpairs[i,] = c('START', lineevents[1])
          i = i + 1
        }
        
        for(j in 2:length(lineevents)-1) # if it's not a terminal, append the event pair
        {
          eventpairs[i,] = c(lineevents[j], lineevents[j+1])
          i = i + 1
        }
        if(terminals) # and append END here
        {
          eventpairs[i,] = c(lineevents[length(lineevents)], 'END')
          i = i + 1
        }
      }
      else if(terminals) # special case where the line is just 1 event, and we need to append START and END
      {
        eventpairs[i,] = c('START', lineevents[1])
        i = i + 1
        eventpairs[i,] = c(lineevents[1], 'END')
      }
      else {
        warning("Single event sequence detected\nRecommend setting terminals to TRUE")
      }
      
    }
    
    # now, if we are collapsing, replace every occurrence of them with the LFE event
    
    if(!is.null(cutoff))
    {
      for(p in 1:nrow(eventpairs)) # for each pair
      {
        for(e in 1:2) # for both of the events in that pair
        {
          if(eventpairs[p,e] %in% lowfreqs) # if that event is low frequency
          {
            eventpairs[p,e] = 'LFE' # replace it with the dummy
          }
        }
      }
    }
    
    return(eventpairs)
  }
  
  getfreqs <- function(seqvec, collapse_arg = collapse) # converts vector of sequence data to transitional frequency matrix
  {  # rows are antecedents and columns are sequiturs
    eventpairs = getpairs(seqvec)
    eventtypes = gettypes(seqvec)
    
    transfreqs = matrix(0, nrow = length(eventtypes), ncol = length(eventtypes), dimnames = list(eventtypes, eventtypes)) # empty trans freq matrix
    
    for(i in 1:nrow(eventpairs))
    {
      #cat("counting transition #", i, "\n")
      if(is.na(eventpairs[i,1]) || is.na(eventpairs[i,2])) # if one of the pair is NA....
      {
        warning("Incomplete event pair at transition #", i, "\nRecommend setting terminals to TRUE")
      }
      else
      {
        transfreqs[eventpairs[i, 1], eventpairs[i, 2]] = transfreqs[eventpairs[i, 1], eventpairs[i, 2]] + 1
      }
    }
    
    # need to remove rows/cols with only zeroes, for chisquare:
    if(length(which(rowSums(transfreqs)==0)) > 0) # but only if there are those rows, else it breaks (this is messy code, thanks R)
    {
      transfreqs = transfreqs[-which(rowSums(transfreqs)==0),]
    }
    if(length(which(colSums(transfreqs)==0)) > 0)
    {
      transfreqs = transfreqs[,-which(colSums(transfreqs)==0)]  
    }
    
    
    return(transfreqs)
  }
  
  getexp <- function(obsfreqs) # gets expected frequency matrix from observed frequencies
  {
    eventtypes = rownames(obsfreqs)
    expfreqs = matrix(0.0, nrow = nrow(obsfreqs), ncol = ncol(obsfreqs), dimnames = list(rownames(obsfreqs), colnames(obsfreqs))) # empty exp freq matrix
    
    for(i in 1:nrow(obsfreqs)) # calculating expected frequencies as ((row total * column total) / grand total)
    {
      for(j in 1:ncol(obsfreqs))
      {
        expfreqs[i,j] = (sum(obsfreqs[i,]) * sum(obsfreqs[,j])) / sum(obsfreqs)
      }
    }
    return(expfreqs)
  }
  

  ##### begin main flow:
  
  observed = getfreqs(seqvec) # observed transitional frequency matrix
  expected = getexp(observed) # expected transitional frequency matrix
  
  getsnr <- function(obsfreqs, expfreqs) # calculate standardised normal residuals based on difference between observed and expected frequencies
  {
    # eventtypes = rownames(obsfreqs)
    snrmatrix = matrix(0, nrow = nrow(obsfreqs), ncol = ncol(obsfreqs), dimnames = list(rownames(obsfreqs), colnames(obsfreqs))) # empty trans freq matrix
    for(i in 1:nrow(snrmatrix))
    {
      for(j in 1:ncol(snrmatrix))
      {
        snr = ((obsfreqs[i,j] - expfreqs[i,j]) / expfreqs[i,j]**0.5) # snr calculated as (cell.obs - cell.exp) / cell.exp**0.5 
        if(!is.nan(snr)) # catch divide by zero errors
        {
          snrmatrix[i,j] = snr
        }
        else
        {
          snrmatrix[i,j] = 0
        }
          
          
      }
    }
    return(snrmatrix)
  }
  
  snrmatrix = getsnr(observed, expected)
  
  ### methods of calculating SNR critical value:

  snrcrit_derivative <- function(snr) # looks for the sharpest turn in the smoothed scree curve
  {
    snrvec = sort(as.vector(snr), decreasing=T)
    snrvec = snrvec[snrvec>0]
    
    secondDerivative <- function(x,i) {return(x[i+1] + x[i-1] - 2 * x[i])}

    snrlen = 1:length(snrvec)
    snrlo = loess(snrvec ~ snrlen, span=0.5)$fitted

    derivatives = numeric()
    
    for (i in 1:length(snrlo)) {derivatives[i] = (secondDerivative(snrlo, i))} # finds max derivative of fitted loess curve, not the original data
    
    plot(snrvec, type="S", ylab="SNR", xlab="Transitions", main="Critical value for standardised residuals (Derivative)")

    lines(snrlo, col="blue", lty=3)
    abline(v=which.max(derivatives), col="red")
    
    return(snrvec[which.max(derivatives)])
  }
  
  snrcrit_projected <- function(snr) # uses vector projection to find the biggest 'crevice' in the scree curve
  {
    snrvec = sort(as.vector(snr), decreasing=T)
    #snrvec_full = snrvec[snrvec>0] # for plotting
    snrvec = snrvec[snrvec>0] # 
    
    origin = c(1, max(snrvec))
    
    endpoint = c(length(snrvec),snrvec[length(snrvec)])
    
    diffs = c(origin[1] - endpoint[1], origin[2] - endpoint[2])
    
    snr_slope = diffs[2]/diffs[1]
    
    #snr_slope = -(max(snrvec)/length(snrvec))
    
    snr_intercept = max(snrvec)-snr_slope
    
    
    distancePointLine <- function(x, y, slope, intercept) { # I borrowed this function from the internet
      # because I don't get geometry 
      # credit to Gregoire Thomas 
      # http://paulbourke.net/geometry/pointlineplane/pointline.r
      ## x, y is the point to test.
      ## slope, intercept is the line to check distance.
      ##
      ## Returns distance from the line.
      ##
      ## Returns 9999 on 0 denominator conditions.
      x1 <- x-10
      x2 <- x+10
      y1 <- x1*slope+intercept
      y2 <- x2*slope+intercept
      distancePointSegment(x,y, x1,y1, x2,y2)
    }
    
    distancePointSegment <- function(px, py, x1, y1, x2, y2) {
      ## px,py is the point to test.
      ## x1,y1,x2,y2 is the line to check distance.
      ##
      ## Returns distance from the line, or if the intersecting point on the line nearest
      ## the point tested is outside the endpoints of the line, the distance to the
      ## nearest endpoint.
      ##
      ## Returns 9999 on 0 denominator conditions.
      lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
      ans <- NULL
      ix <- iy <- 0   # intersecting point
      lineMag <- lineMagnitude(x1, y1, x2, y2)
      if( lineMag < 0.00000001) {
        warning("short segment")
        return(9999)
      }
      u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
      u <- u / (lineMag * lineMag)
      if((u < 0.00001) || (u > 1)) {
        ## closest point does not fall within the line segment, take the shorter distance
        ## to an endpoint
        ix <- lineMagnitude(px, py, x1, y1)
        iy <- lineMagnitude(px, py, x2, y2)
        if(ix > iy)  ans <- iy
        else ans <- ix
      } else {
        ## Intersecting point is on the line, use the formula
        ix <- x1 + u * (x2 - x1)
        iy <- y1 + u * (y2 - y1)
        ans <- lineMagnitude(px, py, ix, iy)
      }
      ans
    }
    
    # back to me:
    
    distances = numeric()
    for(i in 1:length(snrvec))
    {
      distances[i] = distancePointLine(i, snrvec[i], snr_slope, snr_intercept)
    }
    plot(snrvec, type="S", ylab="SNR", xlab="Transitions", main="Critical value for standardised residuals (Vector projected)")


    # fun graph lines:
    abline(snr_intercept, snr_slope, col="blue")
    
    critpoint = c(which.max(distances), snrvec[which.max(distances)])

    b_slope = (origin[2] - critpoint[2]) / (origin[1] - critpoint[1])
    b_intercept = max(snrvec)-b_slope # we end up not drawing this one, but it's the green line
    
    projected_slope = 1/-snr_slope
    projected_intercept = critpoint[2] - (projected_slope*critpoint[1])

    # intersect of snr and projected lines. trust me:
    intersect = c((snr_intercept-projected_intercept)/(projected_slope-snr_slope), 
                  ((projected_slope*snr_intercept)-(snr_slope*projected_intercept)) / (projected_slope-snr_slope))
    
    segments(origin[1], origin[2], critpoint[1], critpoint[2], col="green")
    #segments(intersect[1], intersect[2], critpoint[1], critpoint[2], col="red")
    abline(projected_intercept, projected_slope, col="red")
    
    
    points(origin[1], origin[2], pch=16)
    points(critpoint[1], critpoint[2], pch=16)
    points(endpoint[1], endpoint[2], pch=16)
    points(intersect[1], intersect[2], pch=16)
    
    
    return(snrvec[which.max(distances)])
  }
  
  snrcrit_chisq <- function(chi2, matsize) # sqrt ( chisquare / number of cells)
  {
    critval = (chi2 / matsize)**0.5
    return(critval)
  }
  
  chisq <- function(obsfreqs, expfreqs) # calculate chi-squared statistic of transitional frequencies
  {
    chi2 = 0
    for(i in 1:nrow(obsfreqs))
    {
      for(j in 1:ncol(obsfreqs))
      {
        chi2 = chi2 + (((obsfreqs[i,j] - expfreqs[i,j])**2) / expfreqs[i,j])
      }
    }
    return(chi2)
  }
    
  # calculate critical value depending on run parameters:
  
  if(critval == 'projected')
  {
    transcrit = snrcrit_projected(snrmatrix)
  }   
  else if(critval == 'derivative')
  {
    transcrit = snrcrit_derivative(snrmatrix)
  } 
  else if(critval == 'chisq')
  {
    transcrit = snrcrit_chisq(chisq(observed, expected), length(observed))
    snrvec = sort(as.vector(snrmatrix), decreasing=T)
    snrvec = snrvec[snrvec>0]
    
    plot(snrvec, type="S", ylab="SNR", xlab="Transitions", main="Critical value for standardised residuals (Chisq)")
    abline(v=max(which(snrvec>=transcrit)), col="red")
  }   
  else
  {
    transcrit = critval
    snrvec = sort(as.vector(snrmatrix), decreasing=T)
    snrvec = snrvec[snrvec>0]
    
    plot(snrvec, type="S", ylab="SNR", xlab="Transitions", main="Critical value for standardised residuals (Manual)")
    abline(h=transcrit, col="red")
  }
  
  
  # now identifying the interesting transitions:
  
  notables <- function(snr, obsfreqs, expfreqs, critval, type='high')
  { # returns notably high or low frequency transitions based on SNR matrix and calculated critical value
    pairnum = 1 # iterator
    notable_transitions = data.frame(matrix(ncol=7, nrow=1), stringsAsFactors=F)
    colnames(notable_transitions) = c('Antecedent', 'Sequitur', 'Observed', 'Expected', 'Residual', 'Proportion (Antecedent)', 'Proportion (Sequitur)')
    notable_matrix=observed
    
    for(i in 1:nrow(snr))
    {
      for(j in 1:ncol(snr))
      {
        if(type=='high')
        {
          if((snr[i,j] >= critval) & (observed[i,j] >= min.obs))
          {
            notable_transitions[pairnum,] = c(rownames(snr)[i], colnames(snr)[j], obsfreqs[i,j], round(expfreqs[i,j],3), round(snr[i,j],3), round((obsfreqs[i,j] / sum(obsfreqs[i,])), 3), round((obsfreqs[i,j] / sum(obsfreqs[,j])), 3))
            pairnum = pairnum + 1
          }
          else
          {
            notable_matrix[i,j] = 0
          }
        }
        else if(type=='low')
        {
          if(snr[i,j] <= -critval)
          {
            notable_transitions[pairnum,] = c(rownames(snr)[i], colnames(snr)[j], obsfreqs[i,j], round(expfreqs[i,j],3), round(snr[i,j],3), round((obsfreqs[i,j] / sum(obsfreqs[i,])), 3), round((obsfreqs[i,j] / sum(obsfreqs[,j])), 3))
            pairnum = pairnum + 1
          }
          else
            {
              notable_matrix[i,j] = 0
            }
        }
       }

    }
    #force data types:
    notable_transitions$Observed = as.integer(notable_transitions$Observed)
    notable_transitions$Expected = as.double(notable_transitions$Expected)
    notable_transitions$Residual = as.double(notable_transitions$Residual)
    return(list(notable_transitions, notable_matrix))
  }
  

    

    
    
  hft_obj = notables(snrmatrix, observed, expected, transcrit, type='high')
  lft_obj = notables(snrmatrix, observed, expected, transcrit, type='low')
  
  hft = hft_obj[[1]]
  lft = lft_obj[[1]]
  #hftmatrix = hft_obj[[2]]
  #lftmatrix = lft_obj[[2]]
  
  if(!exists('lowfreqs'))
  {
    lowfreqs = NULL  # we want to output the vector of low frequency events, but if it doesn't exist, we define it as NULL here
  }
  
  eesa_out = list(
      input = seqvec,
      eventsep = eventsep,
      seqsep = seqsep,
      events = countevents(seqvec),
      lfe = lowfreqs,
      observed = observed,
      expected = expected,
      residuals = snrmatrix,
      scree = sort(as.vector(snrmatrix), decreasing=T),
      critical = transcrit,
      hft = hft[order(-hft$Residual),],
      lft = lft[order(-lft$Residual),]
      )
  
  class(eesa_out) = 'eesa'
  
  return(eesa_out)
}

summary.eesa <- function(self) # placeholder
{
  cat('High frequency transitions:\n')
  print(self$hft)
}


eesa.recombine <- function(eesa_obj, critval='default', min.obs='default', junksep='@', hftsep='-')
{
  ordered_hft = eesa_obj$hft
  newdata = eesa_obj$input
  eventsep = eesa_obj$eventsep
  
  if(critval == 'default') {critval = eesa_obj$critical*2} # by default the recombination threshold is twice the critical value for notability
  
  if(min.obs == 'default') {min.obs = 2*mean(eesa_obj$observed[eesa_obj$observed > 0])} # default minimum observations is twice the mean of nonzero observed transitional frequencies. I don't know if this is good or not
  
  for(rownum in 1:nrow(ordered_hft)) # remove the spaces between high frequency events:
  {
    if((ordered_hft$Residual[rownum] > critval) & (ordered_hft$Observed[rownum] > min.obs))
    {
      original = paste(ordered_hft[rownum,1], ordered_hft[rownum,2], sep=eventsep)
      placeholder = paste(junksep, ordered_hft[rownum,1], hftsep, ordered_hft[rownum,2], junksep, sep='') # intermediate step to avoid over-compression
      recombined = paste(ordered_hft[rownum,1], ordered_hft[rownum,2], sep='')
      for(seqnum in 1:length(newdata))
      {
        newdata[seqnum] = gsub(original, placeholder, newdata[seqnum])
      }

    }

  }
  for(seqnum in 1:length(newdata))
  {
    newdata[seqnum] = gsub('@', '', newdata[seqnum])
  }
  return(newdata)
}

eesa.interactive <- function()
{
  cat("Sequence data should be in the form:
      a b c d
      a d d c b
      a b d
      etc.\n")
  cat("Please locate your sequence data file\n")
  dataset = file.choose()
  
  data_terminals = readline("Does your data contain 'start' and 'end' events? > ")
  if(any(grep("[yY]", data_terminals))) {data_terminals = F} else {data_terminals = T} # this is a messy hack, i know
  
  #data_eventsep = readline("Are your events separated by tabs or spaces? > ")
  #if(any(grep("[tT]", data_eventsep))) {data_eventsep = "\t"} else {data_eventsep = " "}
  
  data_min.obs = as.numeric(readline("How many observations should be required for an important transition? "))
  
  cat("Do you want to collapse low frequency events? \n(type 'yes' for default collapsing, or a number to treat that as an observed frequency cutoff)\n")
  data_collapse = readline("> ")
  if(suppressWarnings(!is.na(as.numeric(data_collapse)))) {data_collapse = as.numeric(data_collapse)} else
    if(any(grep("[yY]", data_collapse))) {data_collapse = 'default'} else {data_collapse = NULL}
  
  output = eesa.default(dataset, terminals=data_terminals, min.obs=data_min.obs, collapse=data_collapse)
  View(output$hft)
  
  data_recombine = readline("Do you want to recombine high frequency transitions? > ")
  if(any(grep("[yY]", data_recombine)))
  {
    data_cutoff = readline("At what SNR cutoff? (type 'default' if not sure) > ")
    if(suppressWarnings(!is.na(as.numeric(data_cutoff))))
    {data_cutoff = as.numeric(data_cutoff)} else {data_cutoff = 'default'}
    newdataset = eesa.recombine(output, critval=data_cutoff)
    
    output = eesa.default(newdataset, datatype='v', terminals=data_terminals, min.obs=data_min.obs, collapse=data_collapse, eventsep = data_eventsep)
    View(output$hft)
  }
  
  return(output)
  
}

output = eesa("/home/andrey/Dropbox/Python/arw_data.txt", collapse=NULL, terminals=T, critval='projected')
View(output$hft)
 
#output = eesa.interactive()

