# Exploratory Event Sequence Analysis, by Andrey Barsky
# Version 0.1
# By default, input for this program should be a text file in the format:
# a b c a b c
# a b b c a b b
# etc.
# But different ways of separating events from each other, and sequences from each other,
# can be set with the eventsep and seqsep arguments.

# Any way of specifying events will do - letters, words, numbers or symbols, and however many,
# but should not be the same as either separator argument (for obvious reasons). 
# The @ sign and the dash - are default operators in the recombine function.

# You can also pass different kinds of data as input using the datatype argument, below.

# datatype is an argument specifying the structure of the data you're passing to the eesa function.
# if it is the default 'filetype' or 'f', it treats the input as a file on your machine and reads it as a table
# if it is 'vector' or 'v', it treats it as a character vector of event sequences
# if it is 'string' or 's', it treats it as an atomic character vector separated by linebreaks

# terminals is an optional argument which appends an implicit START and END event to every row of events.
# It is recommended you set it to TRUE, unless your data already has explicit start and end events, 
# or you don't wish to have them for some other reason.

# collapse is an optional argument that allows for collapsing low-frequency events into a single dummy event.
# This can be helpful, since low frequency events are not always interesting but the parts
# of the sequence in which they occur might have some significance.
# By default it is set to NULL, which does no collapsing, but you can set it to 'default'
# to collapse events that occur at least 1 stdev less frequently than the mean, or
# set it to any integer to treat that as the minimum observed frequency cutoff.

# critval specifies the critical value of standardised normal residual
# that marks particular transitions as notable. By default, it is calculated
# as sqrt(chisquare / number of cells), but you can set it to any valid float.

# In the transition matrices, the Y axis lists antecedents. X axis lists consequents.

# Output is an 'eesa' object with the following elements:
# $input: How we parsed the data provided (before collapsing, etc.)
# $alphabet: The unique event codes detected
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


# enter the filepath of your own data here: (I make it try two different paths depending on which computer I'm on)
try((eesa_example = '/home/andrey/Dropbox/Python/arw_data.txt'))
try((eesa_example = '/Users/lpxmac/Dropbox/R/Hayley-Data.txt'))


eesa <- function(x, ...) UseMethod('eesa')

eesa.default <- function(x, datatype='filename', terminals=F, eventsep=' ', seqsep='\n', critval='default', min.obs=1, collapse='default')
{

  ##### different ways of parsing the input:
  
  if(datatype == 'filename' | datatype == 'file' | datatype == 'f')
  {
    seqvec = read.table(x, header=F, sep=seqsep, stringsAsFactors=F)[,1]  
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
    cat("Invalid datatype argument")
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
  }
  else
  {
    cutoff=NULL
  }
  
  if(!is.null(cutoff)) # calculate a vector of low frequency events if we are collapsing
  {
    lowfreqs = lfe(seqvec, cutoff)
  }
  
  
  ##### core functions:
  
  gettypes <- function(seqvec, cutoff_arg = cutoff) # converts vector of sequence data into vector of discrete event types
  {
    alphabet = c()
    for(i in 1:length(seqvec))
    {
      linetypes = strsplit(seqvec[i], eventsep) # break up each line
      alphabet = append(alphabet, linetypes[[1]])
    }
    eventtypes = levels(factor(alphabet)) # get the unique events
    
    if(!is.null(cutoff)) # if we are collapsing, remove the low freq events, and add the LFE event
    {
      eventtypes[!eventtypes %in% lowfreqs] # remove low freq events
      eventtypes = append(eventtypes, 'LFE') # add LFE event
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
  
  
  getpairs <- function(seqvec, cutoff_arg = cutoff) # takes sequence data vector, outputs data frame of event pairs
  {
    eventpairs = data.frame(matrix(ncol=2, nrow=1)) # a growing df with 2 cols
    i = 1
    
    for(line in 1:length(seqvec))
    {
      lineevents = strsplit(seqvec[line], ' ')[[1]]
      
      for(e in lineevents) # fixes errors with doubled separators
      {
        if(e=="")
        {lineevents = lineevents[!lineevents %in% e]}
      }
      
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
  {  # rows are antecedents and columns are subsequents
    eventpairs = getpairs(seqvec)
    eventtypes = gettypes(seqvec)
    #View(eventpairs)
    #View(eventtypes)
    
    transfreqs = matrix(0, nrow = length(eventtypes), ncol = length(eventtypes), dimnames = list(eventtypes, eventtypes)) # empty trans freq matrix
    #View(transfreqs)
    
    for(i in 1:nrow(eventpairs))
    {
      #cat("counting transition #", i, "\n")
      transfreqs[eventpairs[i, 1], eventpairs[i, 2]] = transfreqs[eventpairs[i, 1], eventpairs[i, 2]] + 1
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
  
  # next - chi square and critical values of SNR
  
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
  
  # chi square value of event transitions:

  
  snrcrit <- function(chi2, matsize) # sqrt ( chisquare / number of cells)
  {
    critval = (chi2 / matsize)**0.5
    return(critval)
  }
  
  transchi = chisq(observed, expected)
  if(critval == 'default')
    {
      transcrit = snrcrit(transchi, length(observed))
    }   else
    {
      transcrit = critval
    }
  
  
  # now identifying the interesting transitions:
  
  notables <- function(snr, obsfreqs, expfreqs, critval, type='high')
  { # returns notably high or low frequency transitions based on SNR matrix and calculated critical value
    pairnum = 1 # iterator
    notable_transitions = data.frame(matrix(ncol=5, nrow=1), stringsAsFactors=F)
    colnames(notable_transitions) = c('Antecedent', 'Subsequent', 'Observed', 'Expected', 'Residual')
    notable_matrix=observed
    
    for(i in 1:nrow(snr))
    {
      for(j in 1:ncol(snr))
      {
        if(type=='high')
        {
          if((snr[i,j] > critval) & (observed[i,j] > min.obs))
          {
            notable_transitions[pairnum,] = c(rownames(snr)[i], colnames(snr)[j], obsfreqs[i,j], expfreqs[i,j], snr[i,j])
            pairnum = pairnum + 1
          }
          else
          {
            notable_matrix[i,j] = 0
          }
        }
        else if(type=='low')
        {
          if(snr[i,j] < -critval)
          {
            notable_transitions[pairnum,] = c(rownames(snr)[i], colnames(snr)[j], obsfreqs[i,j], expfreqs[i,j], snr[i,j])
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
  
  eesa_out = list(input = seqvec,
      eventsep = eventsep,
      seqsep = seqsep,
      alphabet = gettypes(seqvec),
      lfe = lowfreqs,
      observed = observed,
      expected = expected,
      residuals = snrmatrix,
      chisq = transchi,
      critical = transcrit,
      hft = hft[order(-hft$Residual),],
      lft = lft[order(-lft$Residual),])
  
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
  
  if(min.obs == 'default') {min.obs = mean(eesa_obj$observed[eesa_obj$observed > 0])} # default minimum observations is the mean of nonzero observed transitional frequencies. I don't know if this is good or not
  
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

# example_vector = read.table(eesa_example, header=F, sep='\n', stringsAsFactors=F)

eesa_output = eesa(eesa_example)

# now you can check out eesa_output and its parameters
