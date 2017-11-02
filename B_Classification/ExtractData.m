function data = ExtractData(filename)

    imported = importdata(filename);

    % PSM1
    n=1;
    for i= 39:41 %  tool position [7]
      data.PSM1.pos(:,n) = imported(:,i);
      n=n+1;
    end

    n=1;
    for i= 42:50 % tool rotations [7]
      data.PSM1.rot(:,n) = imported(:,i);
      n=n+1;
    end

    n=1;
    for i= 51:53 %  linear velocity [7]
      data.PSM1.vel(:,n) = imported(:,i);
      n=n+1;
    end

    n=1;
    for i= 54:56 %  rot. velocity + gripper angle vel. [7]
      data.PSM1.rot_vel(:,n) = imported(:,i);
      n=n+1;
    end
    data.PSM1.grip_vel(:,1) = imported(:,57);
    %% PSM 2

    n=1;
    for i= 58:60 %  tool position [7]
      data.PSM2.pos(:,n) = imported(:,i);
      n=n+1;
    end

    n=1;
    for i= 61:69 % tool rotations [7]
      data.PSM2.rot(:,n) = imported(:,i);
      n=n+1;
    end

    n=1;
    for i= 70:72 %  linear velocity [7]
      data.PSM2.vel(:,n) = imported(:,i);
      n=n+1;
    end

    n=1;
    for i= 73:75 %  rot. velocity + gripper angle vel. [7]
      data.PSM2.rot_vel(:,n) = imported(:,i);
      n=n+1;
    end
    data.PSM2.grip_vel(:,1) = imported(:,76);
end

