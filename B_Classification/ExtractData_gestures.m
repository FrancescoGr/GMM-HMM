function data = ExtractData_gestures(filename)

    imported = importdata(filename);
      data.start = imported(:,1);
      data.stop = imported(:,2);
      data.gestures = imported(:,3);
 
end

