"""
@authors: Juan L. Trincado/Marc Elosua-Bayes
@email: juanluis.trincado@upf.edu/elosua.marc@gmail.com
copy_lims_files.py: this script takes as input the info.txt file from the lims and copies all the files in the output_dir provided, it also takes a metadata file containing the gem_ID mapper between. Marc Elosua adapted the script to incorporate the metadata since he will be runing it on spaceranger.
"""

import os, sys, logging, re, wget
from shutil import copyfile

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

# sample_metadata = "../data/sample_id.txt"
# ids_path = "info.txt"
# output_path = "test2"
def main(ids_path, output_path, sample_metadata):
    try:
        
        logger.info('Starting execution')
        
        path='/scratch/project/production/fastq'
        
        # Create dictionary to save metadata info
        meta_dic = {}
        #1. Load the metadata file (csv) and save it in a dictionary
        with open(sample_metadata, 'r') as f1:
            logger.info('Here b')
            # Skip header
            f1.readline()

            for line1 in f1:
                tokens1 = line1.rstrip().split(',')
                
                # 2nd level index - SI-TT-A1
                ident = tokens1[5]
                meta_dic[ident] = {}
                
                # 3rd fill up the library
                meta_dic[ident]['subproject'] = tokens1[0]
                meta_dic[ident]['gem_id'] = tokens1[1]
                meta_dic[ident]['lib_id'] = tokens1[2]
                meta_dic[ident]['hash'] = tokens1[3]
                meta_dic[ident]['donor_id'] = tokens1[4]
                meta_dic[ident]['slide'] = tokens1[6]
                meta_dic[ident]['area'] = tokens1[7]
                
                # Check if spot layout file exists, if it doesn't, download it!
                # spaceranger requires this files and attempts to download it from
                # the internet. Since jobs cannot access it we need to download it beforehand
                # Link in question is of this format: 
                # https://s3.us-west-2.amazonaws.com/10x.spatial-slides/gpr/V19S23/V19S23-039.gpr
                # Check that it refers to the slide ID
                slide_dir = '../data/spot_layout'
                os.makedirs(slide_dir, 0o755, exist_ok=True);
                fn = '{}.gpr'.format(tokens1[6])
                fn_path = '{}/{}'.format(slide_dir, fn)
                logger.info('Check if spot layout file exists '+fn_path+'...')
                if not os.path.isfile(fn_path):
                    # Split slide ID by -
                    slide_ls = tokens1[6].split('-')
                    url = 'https://s3.us-west-2.amazonaws.com/10x.spatial-slides/gpr/{}/{}-{}.gpr'.format(slide_ls[0],slide_ls[0],slide_ls[1])
                    # Download file
                    logger.info('Downloading '+url+'...')
                    wget.download(url, slide_dir)
                    logger.info('Downloaded '+url+' successfully...')

        
        logger.info('Moving on to writing symbolic links...')
        #2. Load the file with the info of lims files
        with open(ids_path, 'r') as f:
            # Skip header
            f.readline()
            for line in f:
                tokens = line.rstrip().split('\t')
                print(tokens)
                subproject = tokens[1]
                sample_name = re.sub('/', '.', tokens[3])
                logger.info('Processing '+sample_name+'...')
                flowcell = tokens[9]
                lane = tokens[10]
                id = tokens[11]
                status = tokens[12]
                gem_id = meta_dic[id]["gem_id"]
                if(status=='pass'):
                    # fastq_path = output_path+'/'+subproject+'/'+'fastq'
                    fastq_path = "{}/{}/fastq/".format(output_path, subproject)
                    os.makedirs(fastq_path, 0o755, exist_ok=True);
                    #Get the link
                    whole_path_rep1 = path + '/' + flowcell + '/' + lane + '/fastq/' + flowcell + '_' + lane + '_' + id + '_1.fastq.gz'
                    whole_path_rep2 = path + '/' + flowcell + '/' + lane + '/fastq/' + flowcell + '_' + lane + '_' + id + '_2.fastq.gz'
                    logger.info(whole_path_rep1)
                    logger.info(whole_path_rep2)
                    logger.info('Reading symbolic links')
                    link1 = os.readlink(whole_path_rep1)
                    link2 = os.readlink(whole_path_rep2)
                    #Process the link information
                    logger.info('Process the link information')
                    real_path1 = '_'.join(link1.split('/')[-1].split('_')[3:])
                    real_path2 = '_'.join(link2.split('/')[-1].split('_')[3:])
                    # Copy the file into the output directory
                    logger.info('Copying '+output_path+'/'+subproject+'/'+'fastq'+'/'+gem_id+'_'+real_path1+'...')
                    # logger.info('Test '+whole_path_rep1+'...')
                    if not os.path.isfile(output_path+'/'+subproject+'/'+'fastq'+'/'+gem_id+'_'+real_path1):
                      # copyfile(whole_path_rep1, output_path+'/'+sample_name+'_'+real_path1)
                      os.symlink(whole_path_rep1, output_path+'/'+subproject+'/'+'fastq'+'/'+gem_id+'_'+real_path1)
                      
                    logger.info('Copying '+output_path+'/'+subproject+'/'+'fastq'+'/'+gem_id+'_'+real_path2+'...')
                    if not os.path.isfile(output_path+'/'+subproject+'/'+'fastq'+'/'+gem_id+'_'+real_path2):
                      # copyfile(whole_path_rep2, output_path+'/'+sample_name+'_'+real_path2)
                      os.symlink(whole_path_rep2, output_path+'/'+subproject+'/'+'fastq'+'/'+gem_id+'_'+real_path2)

        logger.info('Done. Exiting program.')

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error('Aborting execution')
        sys.exit(1)

# Run function
if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])
