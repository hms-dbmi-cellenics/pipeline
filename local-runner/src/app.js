const AWS = require('aws-sdk');
const fs = require('fs');
const Docker = require('dockerode');
const DockerEvents = require('docker-events');
const { pt } = require('prepend-transform');
const chalk = require('chalk');

// Enable & connect to local Inframock.
AWS.config.update({
  endpoint: 'http://localhost:4566',
  sslEnabled: false,
  s3ForcePathStyle: true,
});

const validPipelineTypes = ['qc', 'gem2s', 'subset', 'copy'];
const isPipelineContainer = (name) => validPipelineTypes.some((keyword) => name.includes(keyword));

const setVarsInTemplate = (template) => {
  const varNames = ['DEBUG_STEP', 'DEBUG_PATH', 'HOST_IP'];
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
    region: 'eu-west-1',
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

  docker.getContainer(id).attach(
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

const attachToNewContainers = (docker, emitter, nameColorMap) => {
  emitter.on('start', (message) => {
    const { id, Actor: { Attributes: { name } } } = message;

    if (!isPipelineContainer(name)) {
      return;
    }
    attachToContainer(docker, id, name, nameColorMap, true);
  });

  const stopDieCallback = (message) => {
    const { Actor: { Attributes: { name } } } = message;

    if (isPipelineContainer(name)) {
      return;
    }

    const coloredName = nameColorMap[name] || name;

    console.log('Container with name', coloredName, 'stopped.');
  };

  emitter.on('stop', stopDieCallback);
  emitter.on('die', stopDieCallback);
  emitter.on('destroy', stopDieCallback);
};

const main = async () => {
  // Initialize stack for Lambdas first.
  try {
    await initStack();
  } catch (e) {
    console.error();
    console.error('Could not initialize stack. Are you sure Inframock is running?');
    console.error();
    process.exit(1);
  }

  // Create hook into Docker API.
  const docker = new Docker();
  const emitter = new DockerEvents({
    docker,
  });

  // Create a color map for easy visualization.
  const nameColorMap = {};

  // Attach to existing running containers and new containers, respectively.
  attachToExistingContainers(docker, nameColorMap);

  attachToNewContainers(docker, emitter, nameColorMap);

  console.log('Waiting for Docker events...');
  emitter.start();
};

main();
