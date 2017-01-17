score_changepoint = function(left,right,smoother=35,eps=.000001){
  left_pgram = get_pgram(left,smoother)
  right_pgram = get_pgram(right,smoother)
  score = mean(left_pgram/(right_pgram+eps)+right_pgram/(left_pgram+eps))
  if (score<0) {
    browser()
  }
  return(score)
}

get_pgram = function(data,smoother=35){
  tmp = fft(data)
  filt = rep(1/smoother,smoother)
  return(filter(Re(tmp*Conj(tmp)),filt,circular=T))
}

score_data = function(data,bandwidth,skip=10,plot_it=TRUE, smoother=35,eps=.000001) {
  scores = numeric(ceiling(length(data)/skip)) #intentionally padded with zeroes to keep bookkeeping easier
  for (i in seq(bandwidth,(length(data)-bandwidth),skip)) {
    if (i%%1000==bandwidth%%1000) {
      print(i/skip)
      if (plot_it){
        plot(scores,type='l')
      }
    }
    left = data[(i-bandwidth+1):i]
    right = data[(i+1):(i+bandwidth)]
    scores[i/skip] = score_changepoint(left,right,smoother,eps)
  }
  if (plot_it){
     plot(scores,type='l')
  }
  return(scores)
}

find_local_peaks = function(scores,window,threshold) {
  if (max(scores)>threshold) {
    peaks = which.max(scores) 
    scores[(peaks-window):(peaks+window)] = 0
  }
  while(max(scores>threshold)) {
    local_max = which.max(scores)
    peaks = c(peaks,local_max)
    scores[(local_max-window):(local_max+window)] = 0
    print(peaks)
  }
  return(peaks)
}
  
rank_frequencies = function(left,right,smoother=11,eps=.0000001) {
  use_length = min(length(left),length(right))
  left = spec.pgram(left[1:use_length],spans=smoother,plot=F)$spec
  right = spec.pgram(right[1:use_length],spans=smoother,plot=F)$spec
  scores = left/(right+eps)+right/(left+eps)
  return(cbind(length(left)/order(scores,decreasing=T),scores[order(scores,decreasing=T)]/mean(scores)))
}
     
grab_region = function(start,finish,skip,trim=.1) {
  start = start*skip
  finish = finish*skip
  offset = (finish-start)*trim
  start = round(start+trim)
  finish = round(finish-trim)
  return(c(start,finish))
} 
