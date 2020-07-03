close all


%%
%path_in = 'E:/AudioSamples/test_scenes/';
%path_render = 'C:/Users/Leo/Documents/renders/';
path_in = '/Users/mccorml1/Documents/AudioSamples/test_scenes/';
path_render = '/Users/mccorml1/Desktop/renders/'; 

audiofiles = {'BAND_shakers_bass_strings_drums.wav', '4SPEAKERS_maleEng_femEng_maleDan_femDan.wav'};
audiofiles_shortnames = {'band', 'speech'};


%% MOVING RECIEVER CASES
order = 3;
for fi = 1:length(audiofiles)
    filepath = [ path_in  audiofiles{fi}];  
    file_shortname = audiofiles_shortnames{fi};
    for si = 4%1:4
        for ri = {'anechoic', 'medium', 'large'} 
            path_out = [path_render file_shortname '_moving_receiver_' num2str(si) 'srcs_' cell2mat(ri) ];
                        if ~exist(path_out, 'dir'), mkdir(path_out); end
            [scene, src_sigs] = A_SETUP_SCENE_moving_receiver(cell2mat(ri), filepath, si, path_out);
            trajectory = B_LISTENER_TRAJECTORY(scene, 1);
            [shmovsig, binmovsig] = C_SYNTHESIZE_MOVING_REF(scene, trajectory, src_sigs, order, 0);
 
            audiowrite([scene.path_out filesep file_shortname '_moving_receiver_' num2str(si) 'srcs_' cell2mat(ri) ...
                '_o' num2str(order) '_ACN_N3D.wav'], shmovsig, scene.fs, 'BitsPerSample', 24); 
            audiowrite([scene.path_out filesep file_shortname '_moving_receiver_' num2str(si) 'srcs_' cell2mat(ri) ...
                '_217kemar.wav'], binmovsig, scene.fs, 'BitsPerSample', 24);
            
            close all

            % save metadata
            path_out = [path_out filesep 'metadata']; %#ok
            if ~exist(path_out, 'dir'), mkdir(path_out); end
            save([path_out filesep 'scene'], 'scene');
            save([path_out filesep 'trajectory'], 'trajectory');
        end
    end
end

return
%% STATIC SOURCE/RECEIVER CASES 
order = 7;
for fi = 1:length(audiofiles)
    filepath = [ path_in  audiofiles{fi}];  
    file_shortname = audiofiles_shortnames{fi};
    for si = 1:4
        for ri = {'anechoic', 'medium', 'large'} 
            clear path_out
            path_out = [path_render file_shortname '_static_' num2str(si) 'srcs_' cell2mat(ri) ];
            if ~exist(path_out, 'dir'), mkdir(path_out); end

            % generate scene save to disk
            [scene, src_sigs] = A_SETUP_SCENE_static_sources(cell2mat(ri), filepath, si, path_out);
            trajectory = B_LISTENER_TRAJECTORY_static_sources(scene, 1);
            [shmovsig, binmovsig] = C_SYNTHESIZE_MOVING_REF_static_sources(scene, trajectory, src_sigs, order, 0);
            audiowrite([scene.path_out filesep file_shortname '_static_' num2str(si) 'srcs_' cell2mat(ri) ...
                '_o' num2str(order) '_ACN_N3D.wav'], shmovsig, scene.fs, 'BitsPerSample', 24);
            audiowrite([scene.path_out filesep file_shortname '_static_' num2str(si) 'srcs_' cell2mat(ri) ...
                '_217kemar.wav'], binmovsig, scene.fs, 'BitsPerSample', 24);
            
            close all
            
            % save metadata
            path_out = [path_out filesep 'metadata']; %#ok
            if ~exist(path_out, 'dir'), mkdir(path_out); end
            save([path_out filesep 'scene'], 'scene');
            save([path_out filesep 'trajectory'], 'trajectory');
        end

    end
end


