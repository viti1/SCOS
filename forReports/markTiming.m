function markTiming(timingFile,startTime)
    T = readtable(timingFile);
    
    ylims = get(gca,'YLim');
    for k=1:size(size(T,1))
        currTime = minutes(T.Time(k) - duration(startTime)); % minutes from recording start
        text(currTime,ylims(2)*0.8,T.Event{k},'Rotation',-90,'Clipping','on','FontSize',12);
    end
    set(gca,'YLim',ylims);
    
    
function nOfMin = difftime(eventTime,startTime)  
    nOfMin  = etime(datevec(datenum(eventTime)),datevec(datenum(startTime)))/60;
    