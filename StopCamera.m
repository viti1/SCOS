   objects = imaqfind;
   for k = 1:numel(objects)
     stop( objects(k) );
     delete(objects(k)); 
   end
msgbox("Turn off Laser");
