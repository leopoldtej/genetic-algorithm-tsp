clear;
clc;
close all;

% Open and read the data file
fileID = fopen("uy734.dat");
rawData = textscan(fileID, "%s%s%s", 50);
fclose(fileID);

% % Video maker (testing)
timestamp = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
% v = VideoWriter(['tsp_animation_', timestamp], 'MPEG-4');
% open(v);

videoWriter = VideoWriter(['bestPathsVideo_', timestamp], 'MPEG-4');
open(videoWriter);

% Mutation Operator (change between 1 to 3 for each mutation type).
mutationOperator = 3;

numCities = length(rawData{1});  % Number of cities
latitudeData = rawData{2};
longitudeData = rawData{3};

% Process the latitude and longitude data
for cityIndex = 1:numCities
    xCoords(cityIndex) = str2num(cell2mat(longitudeData(cityIndex)));
    yCoords(cityIndex) = str2num(cell2mat(latitudeData(cityIndex)));
end

% Normalize the coordinates
normalizedX = (xCoords-min(xCoords))/(max(xCoords)-min(xCoords));
normalizedY = (yCoords-min(yCoords))/(max(yCoords)-min(yCoords));

% Combine the coordinates into one matrix
cityCoordinates = [normalizedX; normalizedY];

% Set up the figure
figure('Units', 'pixels', 'PaperPositionMode','auto');
axis image;
hold on;
plot(normalizedX, normalizedY, 'r*', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
grid;

% Initialize population size and population matrix
populationSize = 1000;  
populationMatrix = zeros(populationSize, numCities);


% Initialize a matrix to store the best paths of specific generations
storedPaths = zeros(7, numCities);
storedGenerations = [1, 50, 100, 500, 1000];

% Populate the population matrix with random routes
for i = 1:populationSize
    populationMatrix(i, :) = randperm(numCities);
end


% Number of generations
numGenerations = 1000;

% Initialize an array to store the shortest distance in each generation
shortestDistances = zeros(numGenerations, 1);

% The Population is represented by the populationMatrix and initialized outside of the given code snippet.
% Each row of the matrix represents an individual in the population.
% Each individual (or chromosome) is a permutation of cities (i.e., a possible route).
% Iterate through each generation
for generation = 1:numGenerations
    % Fitness Function or Evaluation Function
    % Calculate fitness of each individual in the population
    fitnessValues = arrayfun(@(index) calculateDistance(populationMatrix(index, :), cityCoordinates), 1:populationSize);
    
    % Parent Selection Mechanism
    % Perform selection based on fitness values
    [selectedPopulation, ~, selectedIndices] = performSelection(populationMatrix, fitnessValues);
    
    % Recombination (Crossover)
    offspringPopulation = zeros(populationSize, size(populationMatrix, 2));
    for i = 1:2:populationSize
        % Select parents
        parent1 = selectedPopulation(randi([1 populationSize/2]), :);
        parent2 = selectedPopulation(randi([1 populationSize/2]), :);
        
        % Apply Partially Matched Crossover (PMX)
        [child1, child2] = PMX(parent1, parent2);
        
        % Create new generation of offspring
        offspringPopulation(i, :) = child1;
        if i+1 <= populationSize
            offspringPopulation(i+1, :) = child2;
        end
    end

    % Mutation
    mutationRate = .1;  % I found this reduces variability in shortest distance graph (not necessary)
    for i = 1:populationSize
        if rand() < mutationRate
            % Perform mutation
            offspringPopulation(i, :) = performMutation(offspringPopulation(i, :), mutationOperator);
        end
    end

    % The new population becomes the offspring
    populationMatrix = offspringPopulation;
    
    % Fitness Function or Evaluation Function
    % Calculate fitness of each individual in the population after mutation
    fitnessValues = arrayfun(@(index) calculateDistance(populationMatrix(index, :), cityCoordinates), 1:populationSize);
    
    % Save the shortest distance of the current generation
    shortestDistances(generation) = min(fitnessValues);

    % After each generation, check if this is a generation we want to store the best path
    % This is for visualization purposes that Wang wants.
    if any(storedGenerations == generation)
        % Find the best individual and store it
        [~, bestIndex] = min(fitnessValues);
        storedPaths(find(storedGenerations == generation), :) = populationMatrix(bestIndex, :);
    end
end

% 
% for genIndex = 1:numel(storedGenerations)
%     if storedGenerations(genIndex) == generation
%         [~, bestIndex] = min(fitnessValues);
%         storedPaths(genIndex, :) = populationMatrix(bestIndex, :);
%         
%         break;
%     end
% end


% Initialize a new figure for visualization
figure;

% Iterate through each stored generation
for i = 1:length(storedGenerations)
    % Get the best route of the current generation
    bestRoute = storedPaths(i, :);

    % Add the first city to the end of the route to complete the loop
    loopedRoute = [bestRoute, bestRoute(1)];

    % Get the coordinates of the cities in the route
    routeCoordinates = cityCoordinates(:, loopedRoute);

    % Plot the cities
    scatter(routeCoordinates(1, :), routeCoordinates(2, :), 'b.');
    hold on;

    % Plot the route
%     plot(routeCoordinates(1, :), routeCoordinates(2, :), 'r-');
    plot(routeCoordinates(1, :), routeCoordinates(2, :), 'r-', 'LineWidth', 10);


    % Highlight the start/end city
    scatter(routeCoordinates(1, 1), routeCoordinates(2, 1), 200, 'y');

    % Calculate the total distance of the best path
    totalDistance = calculateDistance(bestRoute, cityCoordinates);

    % Set the title with the total distance
    title(['Best Path of Generation ', num2str(storedGenerations(i)), ' - Distance: ', num2str(totalDistance)]);

    % Enable grid
    grid on;

    % Write the frame to the video file
    frame = getframe(gcf);
    for j = 1:60  % Adjust this value to increase or decrease the delay
        writeVideo(videoWriter, frame);
    end

    % Pause for 2 seconds
    pause(2);

    % Clear the figure
    clf;
end


% Close the video file
close(videoWriter);

% Close the figure
close;



% ========= testing visuals ==========
% 
% % Initialize an array to store the final paths
% finpaths = zeros(populationSize, 1);
% 
% % We only want to visualize the last generation, so put the drawing part after the generation loop
% if generation == numGenerations
% 
%     % Set the color map for the routes
%     routeColors = jet(populationSize);
% 
%     % Initialize cell array to store handles of the line segments of each route
%     lineSegmentHandles = cell(1, numCities);
%     
%     % Initialize variable to store the shortest path and its index
%     shortestPath = Inf;
%     shortestPathIndex = 0;
% 
%     % Iterate through each route in the population
%     for routeIndex = 1:populationSize
% 
%         % Get the current route from the population
%         currentRoute = populationMatrix(routeIndex, :);
%     
%         % Create a version of the current route that returns to the start
%         loopedRoute = [currentRoute, currentRoute(1)];
% 
%         % Calculate the total distance of the route
%         totalRouteDistance = calculateDistance(loopedRoute, cityCoordinates);
%         
%         % Store the total distance in finpaths array
%         finpaths(routeIndex) = totalRouteDistance;
%         
%         % Update shortest path if current one is shorter
%         if totalRouteDistance < shortestPath
%             shortestPath = totalRouteDistance;
%             shortestPathIndex = routeIndex;
%         end
% 
%         title(['Distance: ', num2str(totalRouteDistance), '- Mutation: ', mutationOperator]);
%     end
% 
%     % Now plot the shortest path
%     bestRoute = populationMatrix(shortestPathIndex, :);
%     loopedRoute = [bestRoute, bestRoute(1)];
% 
%     % Add mutation operator type to title
%     mutationOperatorType = '';
%     if mutationOperator == 1
%         mutationOperatorType = 'Swap';
%     elseif mutationOperator == 2
%         mutationOperatorType = 'Insert';
%     else
%         mutationOperatorType = 'Scramble';
%     end
% 
%     title(['Best Path Distance: ', num2str(shortestPath), '- Mutation: ', mutationOperatorType]);
% 
%     
% 
%     % Get the coordinates of the start/end city
%     startEndCity = cityCoordinates(:, loopedRoute(1));
% 
%     % Initialize the scatter plot for the start/end city
%     startEndCityScatter = scatter(startEndCity(1), startEndCity(2), 200, 'y');
% 
% 
%     for cityIndex = 1:length(loopedRoute)-1
%         % Plot the segment from the current city to the next one and store the handle
%         line(cityCoordinates(1, loopedRoute(cityIndex:cityIndex+1)), cityCoordinates(2, loopedRoute(cityIndex:cityIndex+1)), 'Color', 'k', 'LineWidth', 2);
%         
% 
%         % Blink the start/end city
%         if mod(cityIndex, 2) == 0
%             set(startEndCityScatter, 'MarkerFaceColor', 'y');
%         else
%             set(startEndCityScatter, 'MarkerFaceColor', 'k');
%         end
% 
%         % Save a video
%         frame = getframe(gcf);  % 'gcf' gets the current figure
%         writeVideo(v, frame);
% 
%         % Write each frame multiple times depending on the desired pause
%         % Here 10 is used assuming your video is at 10 frames per second and you want to pause for 1 second
%         for k = 1:4
%             writeVideo(v, frame);
%         end
% 
%             pause(0.1);
%      end
% 
% end


% close(v);

mutationOperatorType = '';
switch mutationOperator
    case 1
        mutationOperatorType = 'Swap';
    case 2
        mutationOperatorType = 'Insert';
    case 3
        mutationOperatorType = 'Scramble'
end

     
% Define the mutation operator type and the shortest distance in the filename
filename = sprintf('Shortest_Distance_vs_Generation_Mutation_%s_Distance_%.2f.png', mutationOperatorType, shortestDistances(end));

% Plot the shortest distance of each generation
figure;
plot(1:numGenerations, shortestDistances, 'LineWidth', 2);
xlabel('Generation');
ylabel('Shortest Distance');
title(['Shortest Distance vs Generation - Mutation: ', mutationOperatorType]);
grid on;

% Save the plot
saveas(gcf, filename);



% Function to calculate the total distance of a route
function totalDistance = calculateDistance(route, cityCoordinates)
    totalDistance = 0;
    for i = 1:(length(route)-1)
        totalDistance = totalDistance + sqrt((cityCoordinates(1,route(i))-cityCoordinates(1,route(i+1)))^2 + (cityCoordinates(2,route(i))-cityCoordinates(2,route(i+1)))^2);
    end
    totalDistance = totalDistance + sqrt((cityCoordinates(1,route(end))-cityCoordinates(1,route(1)))^2 + (cityCoordinates(2,route(end))-cityCoordinates(2,route(1)))^2);
end
    


% Function to perform selection based on fitness (truncation selection ~ top 50%)
function [newPopulation, newFitnessValues, selectedIndices] = performSelection(population, fitnessValues)
    % Determine the number of individuals to select
    populationSize = size(population, 1);
    numToSelect = floor(populationSize / 2);

    % Sort the fitness values and get the indices in ascending order
    [sortedFitnessValues, sortedIndices] = sort(fitnessValues);
    
    % Select the top 50% individuals
    selectedIndices = sortedIndices(1:numToSelect);
    newPopulation = population(selectedIndices, :);
    newFitnessValues = fitnessValues(selectedIndices);
end


% Function to perform mutation
function individual = performMutation(individual, mutationOperator)   
    switch mutationOperator
        case 1  % Swap
            chrom_size = length(individual);  
            p1 = randi(chrom_size);
            p2 = randi(chrom_size);

            while p1 == p2
                p2 = randi(chrom_size);
            end

            temp = individual(p1);  % use p1 
            individual(p1) = individual(p2);  % use p2 
            individual(p2) = temp;

        case 2  % Insert
            chrom_size = length(individual);
            p1 = randi(chrom_size);
            p2 = randi(chrom_size);
            
            % Make sure p2 is greater than p1
            while p1 >= p2
                p1 = randi(chrom_size);
                p2 = randi(chrom_size);
            end
            
            % Perform the allele shifting
            for i = p2:-1:p1+1
                temp = individual(i);
                individual(i) = individual(i-1);
                individual(i-1) = temp;
            end
    
        case 3  % Scramble
            chrom_size = length(individual);
            scramblePositions = randperm(chrom_size, 2);
            scrambleRange = min(scramblePositions):max(scramblePositions);
            
            % Generate a permutation of the range
            scrambledIndices = randperm(length(scrambleRange));
            
            % Apply the scrambling
            individual(scrambleRange) = individual(scrambleRange(scrambledIndices));
    end
end




% PMX function
function [child1, child2] = PMX(parent1, parent2)
    len = length(parent1);

    % Generate two random crossover points
    crossPoints = sort(randperm(len,2));

    % Copy the genes between crossover points from Parent 1 to child 1
    % and Parent 2 to child 2
    child1 = parent1;
    child2 = parent2;
    
    % Create the mapping between parents
    map1 = parent1(crossPoints(1):crossPoints(2));
    map2 = parent2(crossPoints(1):crossPoints(2));
    
    % Iteratively map the genes not in the crossover section for child1
    for i = [1:(crossPoints(1)-1), (crossPoints(2)+1):len]
        while ismember(child1(i), map1)
            j = find(map1 == child1(i));
            k = find(parent2 == map2(j));
            if child1(k) == map1(j)
                child1(i) = map2(j);
            else
                child1(i) = child1(k);
            end
        end
    end

    % Iteratively map the genes not in the crossover section for child2
    for i = [1:(crossPoints(1)-1), (crossPoints(2)+1):len]
        while ismember(child2(i), map2)
            j = find(map2 == child2(i));
            k = find(parent1 == map1(j));
            if child2(k) == map2(j)
                child2(i) = map1(j);
            else
                child2(i) = child2(k);
            end
        end
    end
end







