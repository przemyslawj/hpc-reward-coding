function [] = playMovies(start,stop, playback_speed, tracking_vid_filename,...
    traces, varargin)

saveVideo = 1;
if saveVideo
    vid_period_ms = 100;
    vid_writer = VideoWriter(['/tmp/movie_' tracking_vid_filename]);
    vid_writer.FrameRate = (1000/200) * 2;
    open(vid_writer);
end

nMov = length(varargin) / 2;
x = ceil(nMov/2);
y = ceil(nMov/x);
f = figure('WindowStyle','normal','Position',[100,200,y*500,x*400]);
frame_ids = ones(nMov, 1);
cvals = [];
for i = 1:nMov
    cvals(i,1) = min(varargin{i * 2 - 1}(:));
    cvals(i,2) = max(varargin{i * 2 - 1}(:));
    stop = min(stop, varargin{i * 2}(end));
end

trace_width = 100;
trace_ahead = 10;
subplot(x,y,nMov+1)
h = plot(double(varargin{4}(1:trace_ahead)) / 1000,...
    createTraceMatrix(traces, 1, trace_ahead));
xlabel('Time (sec)')

for time = start:stop
    if ~ishandle(f)
        break
    end
    for j = 1:nMov
        frame_id = frame_ids(j);
        timestamp = varargin{j * 2}(frame_id);
        if timestamp == time
            subplot(x,y,j)
            movie = varargin{j * 2 - 1};
            imagesc(movie(:,:,frame_id));
            caxis(cvals(j,:));
            colormap(gray)
            
            if j == nMov-1 % ca img movie
                subplot(x,y,nMov+1)
                min_frame = max(1, frame_id - trace_width);
                max_frame = min(size(traces,2), frame_id + trace_ahead);

                A = createTraceMatrix(traces, min_frame, max_frame);
                timestamps = double(varargin{j * 2}(min_frame:max_frame)) / 1000;
                for cell = 1:numel(h)
                    set(h(cell), ...
                        'XData', timestamps,...
                        'YData', A(cell,:))
                end

                drawnow
                
            end
            
            frame_ids(j) = min(size(movie, 3), frame_id + 1);
        end
    end
    pause(1/1000/playback_speed);
    
    if saveVideo && mod(time, vid_period_ms) == 0
        F = getframe(gcf);
        writeVideo(vid_writer, F.cdata)
    end
end

if saveVideo
    close(vid_writer)
end

end

function [A] = createTraceMatrix(traces, start, stop)
    A = traces(:,start:stop);
    for cell = 1:size(A,1)
        A(cell,:) = A(cell,:) + cell;
    end
end
