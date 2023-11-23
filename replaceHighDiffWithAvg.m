function arr = replaceHighDiffWithAvg(arr, threshold1, thresold2)
    % Check if the array has more than 5 elements
    if length(arr) < 5
        error('Array should have at least 5 elements.');
    end

    % Calculate the difference between consecutive elements
    diffs = abs(diff(arr));

    % Find the indices where the difference exceeds the threshold
    highDiffIndices = find(abs(diffs) > threshold1) + 1; % +1 because diff reduces the array size by 1

    % For each of the high difference indices, replace with average of last 5
    for idx = highDiffIndices'
        if idx > 5 
              arr(idx-9:idx+15) = 0; %(arr(idx-5)+arr(idx+5))/2;
        else
            % If there are not enough elements before the index, just use as many as are available
            arr(idx) = mean(arr(1:idx-1));
        end
    end
    diffs = diff(arr);
    highDiffIndices = find((diffs) > thresold2) + 1;
    for idx = highDiffIndices'
        if idx > 5 
              arr(idx-5:idx+10) = 0; %(arr(idx-5)+arr(idx+5))/2;
        else
            % If there are not enough elements before the index, just use as many as are available
            arr(idx) = mean(arr(1:idx-1));
        end
    end

end

