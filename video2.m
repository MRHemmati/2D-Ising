fid = fopen('snapshot.txt');      % open snapshot.txt file; fid is file handle
L = fscanf(fid, '%u', 1);         % find L in 1st line of fid

v = VideoWriter('spinvideo.avi');
open(v);
count = 0;
while (~feof(fid))              % while loop upto the end of fid
  img = ones(L, L);               % generate L*L matrix, which all cells are filled with one
  for i = 1:L
    r = fscanf(fid, '%u', L);     % read one row of image from fid
    if( size(r) == 0)
        break;
    end
    img(i,:) = r;                 % copy r to the i'th image row
  end
  imshow(img);                    % show a frame
  F = getframe(gcf);
  writeVideo(v,F);
  count = count + 1
end
close(v);
fclose(fid);                      % close fid
