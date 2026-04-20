Offset = [15 4 5 ;
          8  5 5 ;
          5  5 15 ] * 10 ;

nFrames = 30;
rec = uint8( round(randn([ size(Offset) , nFrames ])*10  +  repmat(Offset,[1 1 nFrames]) ) ) ;    

WriteAvi('Noise Demo.avi',imresize(rec,30,'nearest'),'Mono8',2);


