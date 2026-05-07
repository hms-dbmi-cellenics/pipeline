const AWS = require('aws-sdk');
const fs = require('fs');
const Docker = require('dockerode');
const { pt } = require('prepend-transform');
const chalk = require('chalk');

// set default values for the AWS region because otherwise setVarsInTemplate function
// will silently fail
process.env.AWS_DEFAULT_REGION = process.env.AWS_DEFAULT_REGION || 'eu-west-1';

// Enable & connect to local Inframock.
AWS.config.update({
  endpoint: 'http://localhost:4566',
  sslEnabled: false,
  s3ForcePathStyle: true,
});

const validPipelineTypes = ['qc', 'gem2s', 'subset', 'obj2s', 'copy'];
const isPipelineContainer = (name) => validPipelineTypes.some((keyword) => name.includes(keyword));

const setVarsInTemplate = (template) => {
  const varNames = ['DEBUG_STEP', 'DEBUG_PATH', 'HOST_IP', 'AWS_DEFAULT_REGION'];
  for (let ii = 0; ii < varNames.length; ii += 1) {
    const value = process.env[varNames[ii]] || '';
    const replace = `__${varNames[ii]}__`;
    const re = new RegExp(replace, 'g');
    // eslint-disable-next-line no-param-reassign
    template = template.replace(re, value);
  }
  return template;
};

const initStack = async () => {
  console.log('Loading CloudFormation for local container launcher...');
  const data = fs.readFileSync('cf-local-container-launcher.yaml', 'utf-8');

  console.log('Creating mock Lambda function on InfraMock...');
  const cf = new AWS.CloudFormation({
    region: process.env.AWS_DEFAULT_REGION,
  });

  const stackName = {
    StackName: 'local-container-launcher',
  };
  try {
    await cf.deleteStack(stackName).promise();
    await cf.waitFor('stackDeleteComplete', stackName).promise();
  } catch (e) {
    console.log('No previous stack found on InfraMock.');
  }

  const { StackId } = await cf.createStack({
    ...stackName,
    Capabilities: ['CAPABILITY_IAM', 'CAPABILITY_NAMED_IAM'],
    TemplateBody: setVarsInTemplate(data),
  }).promise();
  await cf.waitFor('stackCreateComplete', stackName).promise();

  console.log('Stack with ARN', StackId, 'successfully created.');
};

const generateRandomColor = (nameColorMap, name) => {
  // Generate random color for output and easier understanding.
  const randomHexValue = () => Math.floor(Math.random() * 255);
  const coloredName = chalk.rgb(randomHexValue(), randomHexValue(), randomHexValue())(name);

  // eslint-disable-next-line no-param-reassign
  nameColorMap[name] = coloredName;

  return coloredName;
};

const attachToContainer = (docker, id, name, nameColorMap, isNew) => {
  const coloredName = generateRandomColor(nameColorMap, name.replace('/', ''));
  console.log(`Container with name ${coloredName} ${isNew ? 'started' : 'already running'}, attaching log...`);

  const container = docker.getContainer(id);
  container.attach(
    { stream: true, stdout: true, stderr: true },
    (err, stream) => {
      stream.pipe(pt(`${coloredName} | `)).pipe(process.stdout);
    },
  );
};

const attachToExistingContainers = (docker, nameColorMap) => {
  docker.listContainers((err, containers) => {
    containers.forEach((info) => {
      const { Names: names, Id: id } = info;
      const name = names.filter((n) => isPipelineContainer(n))[0];

      if (!name) {
        return;
      }
      attachToContainer(docker, id, name, nameColorMap, false);
    });
  });
};

const attachToNewContainers = (docker, nameColorMap) => {
  docker.getEvents({ filters: { type: ['container'] } }, (err, stream) => {
    let buffer = '';
    
    stream.on('data', (chunk) => {
      buffer += chunk.toString();
      
      const lines = buffer.split('\n');
      buffer = lines[lines.length - 1];
      
      for (let i = 0; i < lines.length - 1; i++) {
        if (!lines[i].trim()) continue;
        
        const message = JSON.parse(lines[i]);
        
        if (message.Type === 'container' && message.Action === 'start') {
          const name = message.Actor.Attributes.name;
          const id = message.Actor.ID;
          
          if (isPipelineContainer(name)) {
            attachToContainer(docker, id, name, nameColorMap, true);
          }
        }
        
        if (message.Type === 'container' && (message.Action === 'stop' || message.Action === 'die' || message.Action === 'destroy')) {
          const name = message.Actor.Attributes.name;
          
          if (isPipelineContainer(name)) {
            const coloredName = nameColorMap[name] || name;
            console.log('Container with name', coloredName, 'stopped.');
          }
        }
      }
    });
  });
};

const main = async () => {
  // Initialize stack for Lambdas first.
  try {
    await initStack();
  } catch (e) {
    console.error();
    console.error('Could not initialize stack. Are you sure Inframock is running?');
    console.error('Actual error:', e.message);
    console.error('Full error:', e);
    console.error();
    process.exit(1);
  }

  // Create hook into Docker API.
  const docker = new Docker();

  // Create a color map for easy visualization.
  const nameColorMap = {};

  // Attach to existing running containers and new containers, respectively.
  attachToExistingContainers(docker, nameColorMap);

  attachToNewContainers(docker, nameColorMap);

  console.log('Waiting for Docker events...');
};

main();
